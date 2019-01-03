
# library imports
library(tidyverse)
library(stringr)
library(limma)
library(psych)
library(Matrix)

# load the raw data: the table starts in Row 5
# there are extra lines at the bottom - need to stop before the end
data_import <- read_tsv("PAW_labeled_grouped_protein_summary_TMT_8.txt", 
                        skip = 4, n_max = 4437, guess_max = 4437)

# the "Filter" column flags contams and decoys
data_all <- filter(data_import, is.na(Filter))

# check number of rows before and after filtering
print("Table lengths before and after filtering:")
cat(nrow(data_import), nrow(data_all))

# make a data frame for each TMT-plex (and drop zeros)
get_data <- function(df, substring) {
    # extracts PAW TMT columns containing a substring
    
    # df: data frame of PAW results
    # substring: a text substring that defines one TMT-plex
    
    # get the columns
    df <- df %>%
      select(starts_with("TotInt")) %>% 
      select(contains(substring))
    
    # drop rows with zeros
    df[apply(df, 1, function(row) all(row !=0 )), ] # return filtered frame
}

# get each experiment and set column names to something easy to parse
A_raw <- get_data(data_all, "_A_")
B_raw <- get_data(data_all, "_B_")
C_raw <- get_data(data_all, "_C_")
D_raw <- get_data(data_all, "_D_")

# automatically find closest channel pairs

SL_norm <- function(df, suffix, print_factors = TRUE) {
    # Normalizes each channel's sum to the average grand total
        # df: data frame of TMT data (one column for each channel)
        # suffix: a text string to append to 1-to-n numbers for column names
        # print_factors: logical to control printing
    
    # compute norm factors and scale columns
    norm_facs <- mean(c(colSums(df))) / colSums(df)
    df_sl  <- sweep(df, 2, norm_facs, FUN = "*")
    
    # simplify column names
    colnames(df_sl) <- str_c(as.character(1:ncol(df)), suffix)
    
    # print the normalization factors for QC check
    if (print_factors == TRUE) {
        cat("\nNormalization Factors: ", suffix, "\n ")
        cat(sprintf("%s - %0.3f\n", colnames(df), norm_facs))
    }
    
    df_sl # return normalized data
}

diff <- function(df, x, y) {
    # computes a standardized difference between two vectors
        # df: data frame that contains x and y
        # x: index of first column
        # y: index of second column
    
    # compute the difference and set the column name
    diff <- 100 * abs(df[x] - df[y]) / rowMeans(df[c(x, y)])
    colnames(diff) <- sprintf("%d_%d", x, y)
    
    diff # return difference vector
}

iqr_grid <- function(df, print_grid = TRUE) {
    # computes an NxN grid of interquartile ranges of diffs    
         # df: a data frame of TMT data
         # print_grid: logical to control printing
    
    # make an empty matrix ncol-by-ncol
    iqrs  <- matrix(0.0, ncol(df), ncol(df))
    
    # populate an upper diagonal "diff" matrix
    for (i in 1:ncol(df)) {
        for (j in (i+1):ncol(df)) {
            if (j > ncol(df)) {
                break
            }
            iqrs[i, j] <- IQR(diff(df, i, j)[[1]])
        }
    }
    
    # print the grid
    if( print_grid == TRUE) {
        cat("\nInterquartile range grid:\n")
        print(round(iqrs, 2))
    }
    
    iqrs # return the iqr grid
}

find_best_pair <- function(df, nmax, suffix, title) {
    # finds channel pairs having the smallest standardized differences    
        # df: data frame with TMT channels
        # nmax: how many channels/pairs to check
        # suffix: suffix to add to channel labels
        # title: title string to use in plots
    
    # first do SL normalization
    df_sl <- SL_norm(df, suffix)
    
    # compute the ird grid
    iqrs <- iqr_grid(df_sl)
    
    # find the nmax number of smallest IQR values
    # from https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
    iqrs_long <- summary(Matrix(iqrs, sparse = TRUE))
    iqrs_ordered <- iqrs_long[order(-iqrs_long[, 3], decreasing = TRUE),]
    candidates <- iqrs_ordered[1, 1] # starting vector of channels
    for (row in 1:nrow(iqrs_long)){
        candidates <- c(candidates, iqrs_ordered[row, 1])
        candidates <- c(candidates, iqrs_ordered[row, 2])
    }
    
    # get the set of indexes and print the best pairs
    for (pair in 1:nmax) {
        top_n <- sort(unique(candidates)[1:nmax])
        i <- iqrs_ordered[pair, 1]
        j <- iqrs_ordered[pair, 2]
        cat(sprintf("\nPair #%d: (%d, %d) => %0.2f", pair, i, j, iqrs[i, j]))
    }
    
    # make the box plots for the top-n pairs
    distributions <- list()
    for (row in 1:nmax){
        distributions[[row]] <- diff(df_sl, iqrs_ordered[row, 1], iqrs_ordered[row, 2])
    }
    boxplot(do.call(cbind, distributions), main = title)
    
    # multi-panel scatter plots for the top-n channels
    cat(pairs.panels(log10(df_sl[top_n]), lm = TRUE, main = title))
}

# nmax is how many top pairs to consider (3 or 4 are recommended values)
top_n <- 4

# find the standards in experiment A
# this is now a one-line function call
find_best_pair(A_raw, top_n, "_A", "Exp A")

# experiment B
find_best_pair(B_raw, top_n, "_B", "Exp B")

# experiment C
find_best_pair(C_raw, top_n, "_C", "Exp C")

# and experiment D
find_best_pair(D_raw, top_n, "_D", "Exp D")

# check the clustering
A_sl <- SL_norm(A_raw, "_A", print_factors = FALSE)
plotMDS(log10(A_sl), main = "Experiment A")

pairs.panels(log10(A_sl), main = "Exp. A - all 10")
