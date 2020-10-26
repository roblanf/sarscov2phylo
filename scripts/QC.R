#!/usr/bin/env Rscript

# usage: Rscript old_tree new_tree threads

args = commandArgs(trailingOnly=TRUE)


# a little QC for the trees
library(ape)
library(phangorn)
library(parallel)
library(ggplot2)

random_nni_distance = function(n, tree){
    # get the distance between tree1 and itself when
    # the second copy has undergone n NNI moves
    
    tree2 = rNNI(tree, moves=n)
    
    return(list(n=n, dist=wRF.dist(tree, tree2)))
}

old = read.tree(args[1])
new = read.tree(args[2])
threads = args[3]

# first we get trees with only the tips they have in common
old_d = drop.tip(old, setdiff(old$tip.label, new$tip.label))
new_d = drop.tip(new, setdiff(new$tip.label, old$tip.label))


obs = wRF.dist(old_d, new_d)

# now we compare this to the same distance calculated at ever increasing numbers of NNI moves
# to get a rough idea of how many random NNI moves this is equivalent to
nnis = seq(from=0, to=10000, length.out=101)

# turn off warnings for this, then turn them back on:
options(warn=-1)
random_nni_dists = mclapply(nnis, random_nni_distance, tree=old_d, mc.cores=threads)
options(warn=0)

# now let's plot the results so we can see what's up
dists = as.data.frame(do.call(rbind, random_nni_dists))
dists$dist = as.numeric(as.character(dists$dist))
dists$n = as.numeric(as.character(dists$n))


thresh = dists$n[max(which(dists$dist<obs))]

p1 = ggplot(dists, aes(x=n, y = dist)) + 
    geom_point() + 
    geom_smooth() + 
    geom_hline(yintercept = obs, size=1, colour='red', alpha = 0.5, linetype='dashed') + 
    ggtitle("Comparison of current tree to tree from previous iteration", 
            subtitle = paste(paste("Observed distance: ", obs, "(red line)"), paste("Equivalent to ~", thresh, ' random NNIs'), sep="\n"))


########
# phylo signal of lineages
# we do this by
# 1. getting the observed tree length of each lineage
# 2. doing the same for 100 random sets of the same number of tips of each lineage
d = read.csv("metadata.csv")

e = subset(d, subset = d$covv_accession_id %in% new$tip.label)
f = e[match(new$tip.label, e$covv_accession_id),]
f$covv_lineage=as.factor(f$covv_lineage)


random_treelen = function(tree, n){
    # choose n random tips and get the tree length
    t = keep.tip(tree, sample(tree$tip.label, n))
    return(sum(t$edge.length))
}

get_signal = function(trait, tree, data){
    # tree is a tree
    # data is a vector IN THE SAME ORDER as the tree tip labels
    # trait is the specific trait in the 'data' vector you want to measure
    
    obs_t = drop.tip(tree, which(data!=trait))
    n = length(obs_t$tip.label)

    obs_l = sum(obs_t$edge.length)
    
    rep_l = sort(replicate(99, random_treelen(tree, n)))
    
    rep_lt = summary(rep_l)
    
    # p value is the position in a ranked list
    p_value = (sum(obs_l>rep_l)+1)/(length(rep_l)+1)

    return(list(trait=as.character(trait), 
                obs = as.numeric(obs_l), 
                p = as.numeric(p_value), 
                n_tips = as.numeric(n),
                rep_Min = as.numeric(rep_lt[1]),
                rep_1stQu = as.numeric(rep_lt[2]),
                rep_Median = as.numeric(rep_lt[3]),
                rep_Mean = as.numeric(rep_lt[4]),
                rep_3rdQu = as.numeric(rep_lt[5]),
                rep_Max = as.numeric(rep_lt[6])))

    
}

clustering = mclapply(levels(f$covv_lineage), get_signal, tree=new, data=f$covv_lineage, mc.cores=threads)
clustering = as.data.frame(do.call(rbind, clustering))
clustering$trait = as.character(clustering$trait)

num = 2:ncol(clustering)
clustering[num] <- sapply(clustering[num],as.numeric)

# only keep the ones where we can reasonably estimate clustering
# so 5 or more members of a lineage
clustering = subset(clustering, n_tips>4)

#### Write out a log file ####
dir.create("./QC")
write.csv(clustering, "./QC/lineage_clustering.csv")
write.csv(dists, "./QC/topological_distances_randomNNI.csv")

pdf("./QC/tree_comparison_plot.pdf")
print(p1)
dev.off()


sink("./QC/QC_report.md", append=TRUE)

cat("# QC Report for the phylogenetic tree\n")
cat("\n")
cat("## Comparison to previous iteration's tree\n")
cat(paste("wRF distance between previous and new iteration's trees:", obs, "\n"))
cat(paste("This is equivalent to a number of random NNIs equal to roughly:", thresh, "\n"))
cat("\n")

cat("## Clustering of Pangolin lineages\n")
cat(paste(sum(clustering$p<0.5), "out of", nrow(clustering), "Pangolin lineages show significant clustering\n"))
cat("\n")
cat("#### Summary of p values\n")

print(summary(clustering$p))
sink()
