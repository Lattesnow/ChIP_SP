# Read the two CSV files
#file1 <- read.delim("*ChIP.xls",sep = "\t")
#file2 <- read.delim("*HiC.xls",sep = "\t")


folder_path <- getwd()
file_list <- list.files(folder_path, pattern = "ChIP\\.xls$", full.names = TRUE)
file_list2 <- list.files(folder_path, pattern = "HiC\\.xls$", full.names = TRUE)
# Initialize an empty list to store the data frames from each file

file1<- read.delim(file_list)
file2<- read.delim(file_list2)


data1 <- file1[rep(1:nrow(file1), each = nrow(file2)), ]

data2 <- file2[rep(1:nrow(file2),nrow(file1)), ]

concatenated_df<- cbind(data1,data2)

# Load the dplyr package
library(dplyr)
# Load the dplyr package

library(dplyr)
filtered_BIN1<- filter(concatenated_df,chr==BIN1_CHR & ((BIN1_START<start & start<BIN1_END)|(BIN1_START<end & BIN1_END>end)))
filtered_BIN2<- filter(concatenated_df,chr==BIN1_CHR & ((BIN2_START<start & start<BIN2_END)|(BIN2_START<end & BIN2_END>end)))

# Create a new bed file with coordinates from HiC matrix
new_df1 <- filtered_BIN1 %>%
  select(chr, BIN2_START, BIN2_END, pileup)
colnames(new_df1)[2] <- "start"
colnames(new_df1)[3] <- "end"

# Create a new data frame with selected columns from data1
new_df2 <- filtered_BIN2 %>%
  select(chr, BIN1_START, BIN1_END, pileup)
colnames(new_df2)[2] <- "start"
colnames(new_df2)[3] <- "end"

final_matrix<- rbind(new_df1,new_df2)

write.table(final_matrix,file = "output_ChIPSP_BedFileA.xls",sep = "\t",row.names = FALSE)
# Print the destination data frame
print(final_matrix)