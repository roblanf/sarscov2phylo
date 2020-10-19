#!/usr/bin/env Rscript

# usage: Rscript prune_terminals.R input_tree output_tree
args = commandArgs(trailingOnly=TRUE)

library(ape)

tree = read.tree(args[1])

drop.terminals <- function(tree, mutations=5){
    
    # iteratively drops terminal branch lengths with >mutation mutations
    # until the tree has no terminal branch lengths longer than this
    # default is 5 mutations
    dropped = data.frame(name = character(), length = numeric())
    
    branchSel = tree$edge[,2]<(Ntip(tree)+1)
    terminal.edgelengths = tree$edge.length[branchSel]
    
    newtree=tree
    while(max(terminal.edgelengths*30000 > mutations)){
        # the true alignment length is just under 30K, so this is equivalent to
        # getting all branches >=mutations mutations
        long.terminals = which(terminal.edgelengths*30000 > mutations)

        newdropped = data.frame(name = newtree$tip.label[long.terminals], length = terminal.edgelengths[long.terminals]*30000)
        dropped = rbind(dropped, newdropped)

        
        # the indices of terminal.edgelengths and tree$tip.label are the same, meaning this works
        newtree = drop.tip(newtree, long.terminals)
        
        branchSel = newtree$edge[,2]<(Ntip(newtree)+1)
        terminal.edgelengths = newtree$edge.length[branchSel]
        
    }
    
    return(list(newtree, dropped))
}

r = drop.terminals(tree)
newtree = r[[1]]
dropped.tips = r[[2]]

if(nrow(dropped.tips)>0){
	dropped.tips$length = paste("terminal branch length of ~", as.integer(dropped.tips[,2]), " mutations", sep="")
}

write.tree(newtree, file=args[2])

write.table(dropped.tips, file="excluded_sequences.tsv", sep="\t", col.names = F, row.names = F, quote=F)
