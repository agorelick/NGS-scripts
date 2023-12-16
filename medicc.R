get_medicc2_tree_and_heatmap <- function(tree_xml_file, cn_profile_file, groups=NULL, title=NA, bp_to_mb=F, vertical=F, ybreaksize=2) { 
    require(ggplot2)
    require(ggtree)
    require(RColorBrewer)
    require(rphyloxml)

    xmltree <- read_phyloxml(tree_xml_file)
    tree <- phyloxml2phylo(xmltree)[[1]]

    ## extract confidence values from xmltree and add to tree
    n_nodes <- length(xmltree[[1]]$node.meta)
    extract_node_info <- function(node, xmltree) {
        node_id <- xmltree[[1]]$edges$id[node]
        l <- unlist(xmltree[[1]]$node.meta[node])
        out <- as.list(xmltree[[1]]$edges[node,])
        toadd <- as.list(l)
        out <- append(out, toadd)
        out
    }
    l <- lapply(1:n_nodes, extract_node_info, xmltree)
    node_dat <- rbindlist(l, fill=T)
    node_dat[isleaf==F, node_label:=name]
    node_dat[!is.na(confidence) & !is.na(node_label),node_label_confidence:=paste0(node_label,' (',confidence,')')]
    node_dat[is.na(confidence) & !is.na(name), node_label_confidence:=name]
    tree$node.label <- node_dat$node_label_confidence[node_dat$isleaf==F]
    tree$tip.label[tree$tip.label=='diploid'] <- 'N'
    tree$node.label <- gsub('internal_','Int_',tree$node.label)

    ## plot the tree+confidences
    if(vertical==T) { 
        p1 <- ggtree(tree,layout='rect',right=F) + geom_nodelab(angle=-90, hjust=-0.2, size=2.5, color='blue') + theme(legend.position='none')
        if(!is.null(groups)) {
            groups[grep('Normal',label), label:='N']
            p1 <- p1 %<+% groups 
            p1 <- p1 + geom_tiplab(aes(color=group), angle=-90, hjust=-0.2) + scale_color_manual(values=group_cols, name='Tissue')
        } else {
            p1 <- p1 + geom_tiplab()
        }
        if(!is.na(title)) p1 <- p1 + ggtitle(label=title) 
        p1 <- p1 + ggplot2::xlim(0, max(p1$data$x)*1.15)
        p1 <- p1 + coord_flip()
        pd <- as.data.table(p1$data)
        pd <- pd[order(y,decreasing=F),]
        pd$x <- max(pd$x) - pd$x
        pd$y[!is.na(pd$label)] <- seq(1,sum(!is.na(pd$label)))
        p1$data <- pd

    } else {
        p1 <- ggtree(tree,layout='rect') + geom_nodelab(angle=0, hjust=0, size=2.5, color='blue') + theme(legend.position='none')
        if(!is.null(groups)) {
            groups[grep('Normal',label), label:='N']
            p1 <- p1 %<+% groups 
            p1 <- p1 + geom_tiplab(aes(color=group)) + scale_color_manual(values=group_cols, name='Tissue')
        } else {
            p1 <- p1 + geom_tiplab()
        }
        if(!is.na(title)) p1 <- p1 + ggtitle(label=title)
        p1 <- p1 + ggplot2::xlim(0, max(p1$data$x)*1.15)
        pd <- as.data.table(p1$data)
        pd <- pd[order(y),]
        pd$y[!is.na(pd$label)] <- seq(1,sum(!is.na(pd$label)))
        p1$data <- pd
    }

    chr <- get_chr_arms()$chr
    scna <- fread(cn_profile_file)
    scna$chrom <- gsub('chr','',scna$chrom)
    if(bp_to_mb) {
        scna[,start:=start/1e6]
        scna[,end:=end/1e6]
    }
    scna[,segid:=paste0(chrom,':',start,'-',end)]
    scna$segid <- factor(scna$segid, levels=unique(scna$segid))
    scna <- merge(scna, chr[,c('chr','global_start'),with=F], by.x='chrom',by.y='chr', all.x=T)
    scna[,global_seg_start:=global_start+start]
    scna[,global_seg_end:=global_start+end]
    scna$sample_id <- gsub('internal_','Int_',scna$sample_id)
    scna[sample_id=='diploid', sample_id:='N']


    ## add in a check for LOH of either allele w.r.t. ancestor node
    tmp <- as.data.table(p1$data)
    tmp <- tmp[,c('parent','node','label'),with=F] 
    tmp$label <- sapply(tmp$label, function(s) strsplit(s,' ')[[1]][1])
    tmp <- merge(tmp[,c('label','node','parent'),with=F], tmp[,c('label','node'),with=F], by.x='parent', by.y='node', all.x=T)
    names(tmp) <- c('parent.node','current.label','current.node','parent.label')
    tmp[,c('current.node','parent.node'):=NULL]
    tmp <- tmp[!is.na(current.label),]
    scna <- merge(scna, tmp, by.x='sample_id', by.y='current.label', all.x=T)
    setnames(scna,'parent.label','parent_id')
    toadd <- scna[,c('segid','sample_id','cn_a','cn_b'),with=F] 
    setnames(toadd,c('cn_a','cn_b'),c('cn_a_anc','cn_b_anc'))
    scna <- merge(scna, toadd, by.x=c('segid','parent_id'), by.y=c('segid','sample_id'), all.x=T)
    scna$is_loh <- F    
    scna[is.na(cn_a_anc), cn_a_anc:=1]
    scna[is.na(cn_b_anc), cn_b_anc:=1]
    scna[(cn_a_anc > 0 & cn_a == 0) | (cn_b_anc > 0 & cn_b == 0), is_loh:=T] ## assume that if the MR ancestor node is not available, the most recent ancestor is the normal

    scna[is_wgd==F & (is_gain==T & is_loss==F), type:='Gain']
    scna[is_wgd==F & (is_gain==F & is_loss==T), type:='Loss']
    scna[is_wgd==F & (is_gain==T & is_loss==T), type:='Both']
    scna[is_wgd==T & (is_gain==F & is_loss==F), type:='WGD']
    scna[is_wgd==T & (is_gain==T & is_loss==F), type:='WGD+Gain']
    scna[is_wgd==T & (is_gain==F & is_loss==T), type:='WGD+Loss']
    scna[is_wgd==T & (is_gain==T & is_loss==T), type:='WGD+Both']

    ## get sample order
    sample_order <- pd$label[!is.na(pd$label)]
    sample_order <- sapply(sample_order, function(s) strsplit(s,' ')[[1]][1], USE.NAMES=F)
    scna$sample_id <- factor(scna$sample_id, levels=sample_order)
    scna$sample_id_int <- as.integer(scna$sample_id)
    scna[is_wgd==T,WGD:='Yes']
    scna[is_wgd==F,WGD:='No']

    cols <- c('#f9d1d1','#d5e3ef','#b7dfb6','#cecece','#ec5e60','#5e97c6','#5eb75c')
    names(cols) <- c('Gain','Loss','Both','WGD','WGD+Gain','WGD+Loss','WGD+Both')

    max_cn <- max(c(scna$cn_a),c(scna$cn_b))
    y_cn=seq(0, max_cn)
    y_breaks <- 0.7*(y_cn / max(y_cn))+0.3/2
    tmp <- data.table(ypos=rep(NA, length(y_breaks)*length(sample_order)))
    sample_num <- rep(1, length(y_breaks))
    global_y_breaks <- y_breaks+rep(0,length(y_breaks))

    for(i in 1:(length(sample_order)-1)) {
        y_cn <- c(y_cn, seq(0, max_cn))
        sample_num <- c(sample_num, rep(i+1, length(y_breaks)))
        global_y_breaks <- c(global_y_breaks, rep(i, length(y_breaks)) + y_breaks)
    }
    tmp$ypos <- global_y_breaks
    tmp$cn <- y_cn
    tmp$sample_num <- sample_num
    if('ypos_a' %in% names(scna)) scna[,ypos_a:=NULL]
    if('ypos_b' %in% names(scna)) scna[,ypos_b:=NULL]

    scna <- merge(scna, tmp, by.x=c('sample_id_int','cn_a'), by.y=c('sample_num','cn'), all.x=T)
    setnames(scna,'ypos','ypos_a')
    scna <- merge(scna, tmp, by.x=c('sample_id_int','cn_b'), by.y=c('sample_num','cn'), all.x=T)
    setnames(scna,'ypos','ypos_b')

    jitter_rows <- which(scna$ypos_a==scna$ypos_b)
    scna[jitter_rows, ypos_a:=ypos_a+0.02]
    scna[jitter_rows, ypos_b:=ypos_b-0.02]

    idvars <- names(scna)
    idvars <- idvars[!idvars %in% c('cn_a','cn_b','ypos_a','ypos_b')]
    scna_a <- scna[,c(idvars,'cn_a','ypos_a'),with=F]
    scna_b <- scna[,c(idvars,'cn_b','ypos_b'),with=F]
    setnames(scna_a,c('cn_a','ypos_a'),c('cn','ypos'))
    scna_a$allele <- 'A'
    setnames(scna_b,c('cn_b','ypos_b'),c('cn','ypos'))
    scna_b$allele <- 'B'
    scna <- rbind(scna_a, scna_b)
    scna$type <- factor(scna$type, levels=c('Gain','Loss','Both','WGD','WGD+Gain','WGD+Loss','WGD+Both'))
    scna[is_loh==T, LOH:='LOH']
    scna[is_loh==F, LOH:='no LOH']
    chr22 <- chr[chr %in% 1:22,]

    allele_cols <- c('black','orange')
    names(allele_cols) <- c('A','B')
    tmp[,ybreaks:=ypos]
    tmp <- tmp[tmp$cn %% ybreaksize == 0]

    if(vertical) {
        p2 <- ggplot(scna) + 
            scale_y_continuous(limits=c(0, max(scna$sample_id_int)), breaks=tmp$ypos, labels=tmp$cn, expand=c(0,0)) + 
            scale_x_continuous(breaks=chr22$global_midpoint, labels=chr22$chr, expand=c(0,0)) + 
            geom_rect(aes(xmin=global_seg_start,xmax=global_seg_end,ymin=sample_id_int-1, ymax=sample_id_int+0, fill=type, linewidth=LOH), color='black') + 
            geom_segment(aes(x=global_seg_start, xend=global_seg_end, y=ypos, yend=ypos, color=allele), linewidth=0.5) +
            geom_hline(yintercept=seq(0,length(sample_order)), color='#bfbfbf', linetype='dashed',size=0.25) +
            geom_vline(xintercept=c(0,chr22$global_end),linewidth=0.25,color='#bfbfbf') +
            scale_color_manual(values=allele_cols, name='Parental copy') + 
            scale_fill_manual(values=cols, na.value='white',name='Change from\nparent node') + 
            scale_linewidth_manual(values=widths, name='LOH since\nparent node?') +
            ang::theme_ang(base_size=9) +
            coord_flip() +
            theme(axis.text.y=element_text(angle=90, size=7,hjust=0.5), axis.line=element_blank(), legend.position='right') + labs(x='Genomic position', y='Copy number')
        p <- plot_grid(p1, p2, ncol=1, align='v', rel_heights=c(2, 4), axis='lr')

    } else {
        p2 <- ggplot(scna) + 
            scale_y_continuous(limits=c(0, max(scna$sample_id_int)), breaks=tmp$ypos, labels=tmp$cn, expand=c(0,0)) + 
            scale_x_continuous(breaks=chr22$global_midpoint, labels=chr22$chr, expand=c(0,0)) + 
            geom_rect(data=scna[is_loh==F], aes(xmin=global_seg_start,xmax=global_seg_end,ymin=sample_id_int-1, ymax=sample_id_int+0, fill=type)) + 
            geom_rect(data=scna[is_loh==T], aes(xmin=global_seg_start,xmax=global_seg_end,ymin=sample_id_int-1, ymax=sample_id_int+0, fill=type), linewidth=0.35, color='black') + 
            scale_color_manual(values=allele_cols, name='Parental copy') + 
            scale_fill_manual(values=cols, na.value='white',name='Change from\nparent node') + 
            geom_segment(aes(x=global_seg_start, xend=global_seg_end, y=ypos, yend=ypos, color=allele), linewidth=0.5) +
            geom_hline(yintercept=seq(0,length(sample_order)), color='#bfbfbf', linetype='dashed',size=0.25) +
            geom_vline(xintercept=c(0,chr22$global_end),linewidth=0.25,color='#bfbfbf') +
            ang::theme_ang(base_size=9) +
            theme(axis.text.y=element_text(angle=90, size=7,hjust=0.5), axis.line=element_blank(), legend.position='right') + labs(x='Genomic position', y='Copy number')
        p <- plot_grid(p1, p2, nrow=1, align='h', rel_widths=c(2, 4), axis='tb')
    }

    p
}




