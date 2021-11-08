#!/usr/local/bin/Rscript
##!/usr/local/R/R-3.6.3/bin/Rscript

print("Load R packages")
library(optparse)
library(DBI)

parse_options = function() {
  option_list = list(
    make_option("--title", help="Plot title (Phenotype, group file)"),
    make_option("--genes_to_plot", help="Gene(s) to plot (may use symbol or Ensembl ID); may give multiple genes separated by comma"),
    make_option("--top_file", help="Gene(s) to plot; to be loaded from 'top results' file"),
    make_option("--top_file_formula", help="Filtering formula for 'top results' file"),
    make_option("--sv_path", help="Path/filename of single variant association results (may use %CHR%)"),
    make_option("--mart_mapping_file", help="File with ENSG mappings", default="mart_export.txt"),
    make_option("--exon_db", help="SQLite3 exon DB, file path; required if no MySQL db given", default=""),
    make_option("--mysql_exon_db", help="MySQL exon DB, database name; required if no SQLite3 DB given", default=""),
    make_option("--mysql_exon_user", help="MySQL exon DB, user name; required if no SQLite3 DB given", default=""),
    make_option("--mysql_exon_password", help="MySQL exon DB, password; required if no SQLite3 DB given", default=""),
    make_option("--mysql_exon_host", help="MySQL exon DB, server host; required if no SQLite3 DB given", default=""),
    make_option("--pdf_output_path", help="Path/filename of PDF output file. May use %SYMBOL%", default="plot_%SYMBOL%.pdf")
  )
  parse_args(OptionParser(option_list=option_list))
}

#################################

prepare_mart_mapping = function(opts) {
  mapping_path = opts$mart_mapping_file
  if (!file.exists(mapping_path)) {
    stop("Need MART mapping file")
  }
  print("====== Read MART mapping table")
  mapping = read.table(mapping_path, h=T, sep="\t")
  colnames(mapping) = c("gene", "symbol", "gene_pos_b38", "gene_chr")
  mapping$gene = as.factor(mapping$gene)
  mapping$symbol = as.factor(mapping$symbol)
  mapping$gene_chr = as.factor(mapping$gene_chr)
  print(paste("Got rows:", nrow(mapping)))
  return(mapping)
}

find_exons_for_gene = function(opts, gene) {
  print("====== Query for exons")
  query = paste("SELECT DISTINCT exon_chrom_start, exon_chrom_end FROM exons WHERE ensembl_gene_id = \"", gene, "\"", sep="")
  if (nchar(opts$mysql_exon_db) > 0) {
    print("Use MySQL backend")
    exon_db = dbConnect(RMySQL::MySQL(), dbname=opts$mysql_exon_db, 
      username=opts$mysql_exon_user, password=opts$mysql_exon_password,
      host=opts$mysql_exon_host)
  } else {
    print("Use SQLite3 backend")
    exon_db = dbConnect(RSQLite::SQLite(), opts$exon_db, synchronous = NULL)
  }
  exons = dbGetQuery(exon_db, query)
  dbDisconnect(exon_db)
  print(paste("Got ", nrow(exons), " exon positions from ",
    min(exons$exon_chrom_start), " to ", max(exons$exon_chrom_end),
    " for ", gene, ".", sep = ""))
  return(exons)
}

find_strand = function(opts, gene) {
  print("====== Strand")
  query = paste("SELECT DISTINCT strand FROM exons WHERE ensembl_gene_id = \"", gene, "\"", sep="")
  if (nchar(opts$mysql_exon_db) > 0) {
    print("Use MySQL backend")
    exon_db = dbConnect(RMySQL::MySQL(), dbname=opts$mysql_exon_db, 
      username=opts$mysql_exon_user, password=opts$mysql_exon_password,
      host=opts$mysql_exon_host)
  } else {
    print("Use SQLite3 backend")
    exon_db = dbConnect(RSQLite::SQLite(), opts$exon_db, synchronous = NULL)
  }
  strand = dbGetQuery(exon_db, query)
  dbDisconnect(exon_db)
  print(paste("Got ", strand, " for ", gene, ".", sep = ""))
  return(strand)
}

find_mapping_for_gene = function(mapping, gene_to_plot) {
  print(paste("Lookup gene:", gene_to_plot))
  mapping2 = mapping[mapping$symbol == gene_to_plot,]
  if (nrow(mapping2) == 0) {
    print("Didn't find symbol, try Ensembl ID")
    mapping2 = mapping[mapping$gene == gene_to_plot,]
  }
  mapping = mapping2
  if (nrow(mapping) != 1) {
    print(paste("Given gene not found/not unique. Got rows:", nrow(mapping)))
  }
  print(mapping)
  return(mapping)  
}

load_sv_for_gene = function(opts, mapping) {
  print("====== Read single-variant file")
  sv_path = opts$sv_path
  if (nchar(sv_path) == 0) {
    stop("Need path to SV association results")
  }
  chr = mapping$gene_chr[1]
  gene = mapping$gene[1]
  symbol = mapping$symbol[1]
  sv_fn = gsub("%CHR%", chr, sv_path)
  print(paste("Using Ensembl gene:", gene))
  print(paste("Symbol:", symbol))
  print(paste("Filename:", sv_fn))
  print("Subset and remove missing betas")
  d = read.table(sv_fn, h=T)
  d$gene = as.factor(d$gene)
  print(paste("All rows:", nrow(d)))
  e = d[as.character(d$gene) == as.character(gene),]
  print(paste("Rows for gene:", nrow(e)))
  e_nonmiss = e[!is.na(e$beta),]
  print(paste("Without missing betas:", nrow(e_nonmiss)))
  e_nonmiss = merge(e_nonmiss, mapping)
  print(paste("Merged to mapping:", nrow(e_nonmiss)))
  print("Determine variant positions")
  for (i in 1:nrow(e_nonmiss)) {
    variant = as.character(e_nonmiss[i, "Name"])
    items = unlist(strsplit(variant, ":", fixed=T)[[1]])
    e_nonmiss[i, "chr"] = items[1]
    e_nonmiss[i, "pos"] = items[2]
    e_nonmiss[i, "all1"] = items[3]
    e_nonmiss[i, "all2"] = items[4]
  }
  e_nonmiss$chr = as.factor(e_nonmiss$chr)
  e_nonmiss$pos = as.numeric(e_nonmiss$pos)
  e_nonmiss$all1 = as.factor(e_nonmiss$all1)
  e_nonmiss$all2 = as.factor(e_nonmiss$all2)
  e_nonmiss$log_p = -log10(e_nonmiss$p)
  print(summary(e_nonmiss))
  return(e_nonmiss)
}

exon_has_variants = function(exon, variants) {
  variant_count_in_exon = length(which(variants$pos >= exon$exon_chrom_start[1] &
                              variants$pos < exon$exon_chrom_end[1]))
  return(variant_count_in_exon > 0)
}

find_first_exon_with_variants = function(exons, variants) {
  i = 1
  while (i <= nrow(exons)) {
    if (exon_has_variants(exons[i,], variants)) {
      return(i)
    }
    i = i + 1
  }
  return(-1)
}

find_last_exon_with_variants = function(exons, variants) {
  i = nrow(exons)
  while (i > 0) {
    if (exon_has_variants(exons[i,], variants)) {
       return(i)
    }
    i = i - 1
  }
  return(-1)
}

trim_exons_without_variants = function(variants, exons) {
  exons = exons[order(exons$exon_chrom_start, exons$exon_chrom_end),]
  last_exon_with_variants = find_last_exon_with_variants(exons, variants)
  first_exon_with_variants = find_first_exon_with_variants(exons, variants)
  if (first_exon_with_variants == -1 || last_exon_with_variants == -1) {
    #print("No exons overlap variants - plot without exons")
    #return(data.frame())
    print("No exons overlap variants")
    return(exons)
  } else {
    print(paste("Trimming exons at the beginning and the end; start with index ", first_exon_with_variants, 
                " and end with index ", last_exon_with_variants, sep=""))
    #return(exons[first_exon_with_variants:last_exon_with_variants,])
    return(exons)
  }
}

find_genes_to_plot_from_top_file = function(opts) {
  top = read.table(opts$top_file, h = T)
  print(paste("Loaded ", nrow(top), " genes from file: ", opts$top_file, sep = ""))
  expr = parse(text = opts$top_file_formula)
  print(paste("Filter top results using expression: ", expr, sep=""))
  top2 = subset(top, eval(expr, envir = top))
  print(paste("Got ", nrow(top2), " filtered results.", sep=""))
  return(top2$gene)
}

do_plot = function(opts, e_nonmiss, exons, strand) {
  print("====== Plot")
  
  strand.dir = strand

  symbol = e_nonmiss[1,"symbol"]
  chrom  = e_nonmiss[1,"chr"]
  out_fn = opts$pdf_output_path
  out_fn = gsub("%SYMBOL%", symbol, out_fn)
  print(paste("Write plot to:", out_fn))
  pdf(out_fn)
  
  x = e_nonmiss$pos/1000
  y = e_nonmiss$beta
  
  x.exon.sta = exons$exon_chrom_start/1000
  x.exon.end = exons$exon_chrom_end/1000

  circle_size = sqrt(e_nonmiss$log_p/pi) # surface-proportional

  maxMinusLog10p = max(e_nonmiss$log_p)
  minMinusLog10p = min(e_nonmiss$log_p)

  minx = min(c(floor(x), floor(exons$exon_chrom_start/1000)))
  maxx = max(c(ceiling(x), ceiling(exons$exon_chrom_end/1000)))

  minx.trim = min(floor(x)) - 0.5
  maxx.trim = max(ceiling(x)) + 0.5

  diffRange.trim = maxx.trim - minx.trim

  miny = min(c(floor(y - e_nonmiss$se), -1)) - 1
  maxy = max(c(ceiling(y + e_nonmiss$se), 1))

  maxEffect = max(e_nonmiss$beta)
  diffRange = maxy - miny
  rectH = 0.01 * abs(diffRange) #0.05 * maxEffect

  exons.trim <- exons
  trim.left  <- 0
  trim.right <- 0
  for (i in 1:length(exons.trim[,1])) {
       if (minx.trim > exons.trim[i,"exon_chrom_start"]/1000 & minx.trim <= exons.trim[i,"exon_chrom_end"]/1000) {
             exons.trim[i,"exon_chrom_start"] <- minx.trim*1000
       }
       if (minx.trim > exons.trim[i,"exon_chrom_end"]/1000) {
             exons.trim[i,"exon_chrom_start"] <- NA
             exons.trim[i,"exon_chrom_end"] <- NA
       }
       if (!is.na(exons.trim[i,"exon_chrom_start"]) & !is.na(exons.trim[i,"exon_chrom_end"]) &
           maxx.trim > exons.trim[i,"exon_chrom_start"]/1000 & maxx.trim <= exons.trim[i,"exon_chrom_end"]/1000) {
             exons.trim[i,"exon_chrom_end"] <- maxx.trim*1000
       }
       if (!is.na(exons.trim[i,"exon_chrom_start"]) & maxx.trim < exons.trim[i,"exon_chrom_start"]/1000) {
             exons.trim[i,"exon_chrom_start"] <- NA
             exons.trim[i,"exon_chrom_end"] <- NA
             trim.right <- 1
       }
  }
  if (minx < minx.trim) trim.left <- 1
  if (maxx > maxx.trim) trim.right <- 1
 
  exons.trim <- exons.trim[!is.na(exons.trim[,"exon_chrom_start"]) & !is.na(exons.trim[,"exon_chrom_start"]),]

  ##
  
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)

  # plot w/o trimming
  
  symbols(x = x, y = y, fg = "white", 
          circles = circle_size,
          inches = 1/8,
          ann = F, 
          xlim = c(minx, maxx),
          ylim = c(miny, maxy),
          bty = "n")

  lines(x = c(minx, maxx), y = c(0, 0), col = gray(level = 0.3), lty = 2)
 
  lines(x = c(min(exons$exon_chrom_start)/1000, 
              max(exons$exon_chrom_end)/1000),
        y = c(miny+rectH/2,miny+rectH/2), col = "black")

  if (nrow(exons) > 0) {
    rectHeight = rectH
    rect(border = "black", col = "black",
         xleft = exons$exon_chrom_start / 1000, xright = exons$exon_chrom_end / 1000,
         ytop = miny+rectHeight, ybottom = miny-rectHeight/2)
  }

  a.x <- minx+(maxx-minx)/5 #min(exons$exon_chrom_start/1000)
  text(x = a.x, y = miny+2*rectHeight, adj = c(0,0), cex = 0.8, 
       labels = paste(symbol), font=3)   
 
  a0 <- minx+(maxx-minx)/2
  a1 <- minx+1.3*(maxx-minx)/2
  b0 <- miny+3*rectHeight
  if (strand.dir ==  1) {
    arrows(x0 = a0, x1 = a1, y0 = b0, y1 = b0, length=0.05, col = gray(level=0.3), lty = 1)
  } else {
    arrows(x0 = a1, x1 = a0, y0 = b0, y1 = b0, length=0.05, col = gray(level=0.3), lty = 1)
  }

  arrows(x, y - e_nonmiss$se, x, y + e_nonmiss$se, length=0.05, angle=90, code=3, col="gray")

  symbols(x = x, y = y, 
          circles = circle_size,
          inches = 1/8,
          ann = F, 
          bg = "steelblue2", 
          fg = "black",
          add = T,
          bty = "n")
  
  title(main = paste(opts$title),
        xlab = paste("chr ",chrom,", position [kb]",sep=""), 
        ylab = "Effect Size",
        sub = paste("N = ", nrow(e_nonmiss), " SNVs", sep=""))
  mtext(symbol, side = 3, line = 0, font = 4, cex = 1.2)
 
  legend_count = 4
  legend_entries = (max(circle_size) / legend_count) * 1:legend_count
  legend_labels = round((legend_entries ^ 2) * pi, 1)
  
  # cex=3.3 makes a circle of 1/8 inch size
  legend_entries = (3.3 / legend_count) * 1:legend_count

  legend(x = "topright",
         ncol = 1,
         pch = 21, 
         pt.cex = legend_entries,
         legend = legend_labels,
         pt.bg = "steelblue2",
         bg = "transparent",
         inset=c(-0.2,0),
         y.intersp = 1.5,
         title = as.expression(bquote(paste('-log'['10']*'(p)')))
  )

  # plot w trimming
  
  symbols(x = x, y = y, fg = "white",
          circles = circle_size,
          inches = 1/8,
          ann = F, 
          xlim = c(minx.trim, maxx.trim),
          ylim = c(miny, maxy),
          bty = "n")

  lines(x = c(minx.trim, maxx.trim), y = c(0, 0), col = gray(level = 0.3), lty = 2)

  lines(x = c(min(exons.trim$exon_chrom_start)/1000, 
              max(exons.trim$exon_chrom_end)/1000),
        y = c(miny+rectH/2,miny+rectH/2), col = "black")

  if (nrow(exons.trim) > 0) {
    rectHeight = rectH
    rect(border = "black", col = "black",
         xleft = exons.trim$exon_chrom_start / 1000, xright = exons.trim$exon_chrom_end / 1000,
         ytop = miny+rectHeight, ybottom = miny-rectHeight/2)
  }
  
  a.x.trim <- minx.trim+(maxx.trim-minx.trim)/5 #min(exons.trim$exon_chrom_start/1000)
  text(x = a.x.trim, y = miny+2*rectHeight, adj = c(0,0), cex = 0.8, 
       labels = paste(symbol), font=3)

  a0.trim <- minx.trim+(maxx.trim-minx.trim)/2
  a1.trim <- minx.trim+1.3*(maxx.trim-minx.trim)/2
  b0.trim <- miny+3*rectHeight
  if (strand.dir ==  1) {
    arrows(x0 = a0.trim, x1 = a1.trim, y0 = b0.trim, y1 = b0.trim, length=0.05, col = gray(level=0.3), lty = 1)
  } else {
    arrows(x0 = a1.trim, x1 = a0.trim, y0 = b0.trim, y1 = b0.trim, length=0.05, col = gray(level=0.3), lty = 1)
  }

  arrows(x, y - e_nonmiss$se, x, y + e_nonmiss$se, length=0.05, angle=90, code=3, col="gray")

  symbols(x = x, y = y, 
          circles = circle_size,
          inches = 1/8,
          ann = F, 
          bg = "steelblue2", 
          fg = "black",
          add = T,
          bty = "n")

  add.line <- 0.03 * abs(diffRange.trim)

  if (trim.left == 1) {
      lines(x = c(min(exons.trim$exon_chrom_start)/1000-add.line, 
                  min(exons.trim$exon_chrom_start)/1000),
            y = c(miny+rectH/2,miny+rectH/2), lty = 3, col = gray(level = 0.3))
  }

  if (trim.right == 1) {
      lines(x = c(max(exons.trim$exon_chrom_end)/1000, 
                  max(exons.trim$exon_chrom_end)/1000+add.line),
            y = c(miny+rectH/2,miny+rectH/2), lty = 3, col = gray(level = 0.3))
  }

  title(main = paste(opts$title),
        xlab = paste("chr ",chrom,", position [kb]",sep=""), 
        ylab = "Effect Size",
        sub = paste("N = ", nrow(e_nonmiss), " SNVs", sep=""))

  mtext(symbol, side = 3, line = 0, font = 4, cex = 1.2)

  legend_count = 4
  legend_entries = (max(circle_size) / legend_count) * 1:legend_count
  legend_labels = round((legend_entries ^ 2) * pi, 1)
  
  # cex=3.3 makes a circle of 1/8 inch size
  legend_entries = (3.3 / legend_count) * 1:legend_count

  legend(x = "topright",
         ncol = 1,
         pch = 21, 
         pt.cex = legend_entries,
         legend = legend_labels,
         pt.bg = "steelblue2",
         bg = "transparent",
         inset=c(-0.2,0),
         y.intersp = 1.5,
         title = as.expression(bquote(paste('-log'['10']*'(p)')))
  )

  dev.off()
}

load_and_plot = function(opts, mapping, gene_to_plot) {
  mapping = find_mapping_for_gene(mapping, gene_to_plot)
  if (nrow(mapping) > 0) {
    e_nonmiss = load_sv_for_gene(opts, mapping)
    exons = find_exons_for_gene(opts, mapping$gene[1])
    exons = trim_exons_without_variants(e_nonmiss, exons)
    strand = find_strand(opts, mapping$gene[1])
    do_plot(opts, e_nonmiss, exons, strand)
  }
}

plot_genes = function(opts) {
  mapping = prepare_mart_mapping(opts)
  genes_to_plot = c()
  if (length(opts$genes_to_plot) > 0) {
    genes_to_plot = trimws(unlist(strsplit(opts$genes_to_plot, ",", fixed=T)[[1]]))
  }
  if (length(genes_to_plot) == 0) {
    genes_to_plot = find_genes_to_plot_from_top_file(opts)
  }
  if (length(genes_to_plot) == 0) {
    stop("Need genes to plot.")
  }
  for (gene_to_plot in genes_to_plot) {
    load_and_plot(opts, mapping, gene_to_plot)
  }
  print("Finished")
}

#################################

main = function() {
  opts = parse_options()
  plot_genes(opts)
}

main()
