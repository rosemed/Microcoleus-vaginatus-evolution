merge_intervals <- function(df) {
  df <- df[order(df$Beg), ]
  merged <- list()
  current_interval <- as.numeric(df[1, 2:3])
  for (i in 2:nrow(df)) {
    interval <- as.numeric(df[i, 2:3])
    if (interval[1] <= current_interval[2]) {
      current_interval[2] <- max(current_interval[2], interval[2])
    } else {
      merged <- c(merged, list(current_interval))
      current_interval <- interval
    }
  }
  merged <- c(merged, list(current_interval))
  result <- do.call(rbind, merged)
  colnames(result) <- c("Beg", "End")
  return(result)
}
merge_and_count_overlap <- function(df) {
  if (nrow(df)<0){
  print(paste0("no recombination overlap in G",i))}else{
  starts <- data.frame(pos = df$Beg, type = 1)
  ends <- data.frame(pos = df$End, type = -1)
  all_points <- rbind(starts, ends)
  all_points <- all_points[order(all_points$pos), ]
  
  current_overlap <- 0
  overlap_counts <- list()
  current_start <- NULL
  
  for (i in 1:nrow(all_points)) {
    point <- all_points[i, ]
    if (is.null(current_start)) {
      current_start <- point$pos
    }
    current_overlap <- current_overlap + point$type
    if (i < nrow(all_points) && point$pos != all_points[i + 1, "pos"]) {
      overlap_counts[[length(overlap_counts) + 1]] <- c(current_overlap,current_start, point$pos)
      current_start <- all_points[i + 1, "pos"]
    }
  }
  
  result <- do.call(rbind, overlap_counts)
  colnames(result) <- c("overlap_count","Beg", "End")
  return(result)
}
}
get_og_combinations <- function(intervals, og_data) {
  og_combinations <- list()
  
  for (i in 1:nrow(intervals)) {
    current_interval_begin <- intervals$Beg[i]
    current_interval_end <- intervals$End[i]
    
    overlapping_rows <- og_data[
      which((og_data$Beg <= current_interval_begin) & (og_data$End >= current_interval_begin)):which((og_data$Beg <= current_interval_end) & (og_data$End >= current_interval_end)),
    ]
    
    overlapping_ogs <- overlapping_rows$OG
    
    og_combination <- paste(overlapping_ogs, collapse = ",")
    
    og_combinations[[i]] <- og_combination
  }
  
  intervals$OG_combination <- unlist(og_combinations)
  
  return(intervals)
}

loc=read.table("OG_location.txt",header=T)
sequence_length =2284707 #length of core single-copy gene concatenated sequence
result=NULL
for (i in c(2,4,5,7,8)){
file=paste0('within_G',i,'.importation_status.txt')
outfile1=paste0('G',i,'_recombination_interval.txt')
outfile2=paste0('G',i,'_recombination_overlaploc.txt')
data=read.table(file,header=T)
merged_data <- as.data.frame(merge_intervals(data))
intervals <- get_og_combinations(merged_data, loc)
total_covered_length <- sum(intervals$End - intervals$Beg + 1)
interval_midpoints <- (intervals$Beg + intervals$End) / 2
average_distance <- mean(diff(interval_midpoints))
coverage <- rep(0, sequence_length)
for (row in 1:nrow(intervals)) {
  coverage[intervals$Beg[row]:intervals$End[row]] <- 1
}
coverage_freq <- table(coverage) / sequence_length
entropy <- -sum(coverage_freq * log2(coverage_freq), na.rm = TRUE)
result=rbind(result,data.frame(phylogroup=paste0('G',i),coverage_length=total_covered_length,average_distance=average_distance,entropy=entropy))
merged_overlap <- as.data.frame(merge_intervals(subset(as.data.frame(merge_and_count_overlap(data)),overlap_count>length(unique(data$Node))/3)))
overlap_og <- get_og_combinations(merged_overlap, loc)
write.table(intervals,outfile1,row.names=F,sep="\t",quote=F)
write.table(overlap_og,outfile2,row.names=F,sep="\t",quote=F)
}
write.table(result,'recombination_distribution.txt',row.names=F,sep="\t",quote=F)