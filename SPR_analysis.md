# How many SPR moves do I need in fasttree?

This is a question I tried to address before, and settled on a (very) conservatively high number of 10. Since then, I've been updating the global tree every two days, by adding seuqences with Maximum Parsimony and updating them the ME and ML in fasttree. By far the slowest step are the ME SPR moves, and I currently do 10 rounds of range 1000 (i.e. the range is the whole tree).

In this document I'll just take a quick look at how many changes ever occur in each of those 10 rounds of SPR moves over about 10 update cycles. To do that, I'll use the logs I recorded from fasttree to just measure Robinson-Foulds distance between the starting tree and each of the 10 ME-SPR trees from each update. I started loggin the fasttree output in the 14-08-20 release. So I'll start with that one here too.

Below is some clunky R code to get the information I want. And here's the result.

![SPR iterations plot](https://github.com/roblanf/sarscov2phylo/blob/master/spr_iterations.png)

There are a LOT of limitations here. RF distances are confounded for many reasons here, but they're the only thing fast enough to practically calculate. Still, the results show clearly that the amount that the tree changes tails off a lot with each iteration. Almost half of all the changes come in the first iteration. 

My strong suspicion is that most of the changes at all stages are small changes near the tips of the trees (see also the previous analyses of that issue on this repo). And combined with this (highly flawed but still informative) analysis, I think it's sensible to drop the number of SPR rounds to 2.



```
library(readr)
library(phangorn)
library(parallel)
library(ggplot2)

# here's a function to get the distance I want from a fastree log file
get_rfdists = function(filename){

    d = read_lines(filename, skip = 7, n_max = 11)
    e = read_delim(d, delim="\t", col_names = FALSE)
    e$RF=NA
    e$RFn=NA
    
    rownum=2:11
    for(i in rownum){
        
        t1 = read.tree(text=e$X2[i-1])
        t2 = read.tree(text=e$X2[i])
        
        e$RF[i] = RF.dist(t1, t2, normalize = FALSE)
        e$RFn[i] = RF.dist(t1, t2, normalize = TRUE)
        
    }

    e$name = basename(dirname(filename))
    
    r1 = data.frame(name = e$name[2:11], dist = e$RF[2:11], dist.type = "rf.dist")
    r1$spr.iteration = 1:10

    r2 = data.frame(name = e$name[2:11], dist = e$RFn[2:11], dist.type = "rf.dist.norm")
    r2$spr.iteration = 1:10
    
    r = rbind(r1, r2)
            
    return(r)
    
}

# here are the log files
files = c("global_14-08-20/fasttree.log",
          "global_16-8-20/fasttree.log",
          "global_18-08-20/fasttree.log",
          "global_20-08-20/fasttree.log",
          "global_22-08-20/fasttree.log",
          "global_24-08-20/fasttree.log",
          "global_26-08-20/fasttree.log",
          "global_28-08-20/fasttree.log",
          "global_30-08-20/fasttree.log",
          "global_30-08-20/fasttree.log",
          "global_02-09-20/fasttree.log",
          "global_04-09-20/fasttree.log",
          "global_06-09-20/fasttree.log",
          "global_08-09-20/fasttree.log",
          "global_10-09-20/fasttree.log")

# run the function and bind the results
r = mclapply(files, get_rfdists, mc.cores = 8)
results = do.call(rbind, r)

# plot it out
ggplot(results, aes(x=spr.iteration, y=dist, colour=name)) + geom_point() + geom_line() + facet_wrap(facets = c("dist.type"), ncol = 1, scales = "free_y")

```