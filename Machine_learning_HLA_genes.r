
#Setting the working directory 
setwd("~/Documents/UOFG/software_tools_class")

#Installing,updating ,and loading libraries 
BiocManager::install(version = "3.15")
BiocManager::install("Biostrings",force = TRUE)
install.packages('bioseq')
install.packages('qpcR')
install.packages("rgl")
install.packages("glmnet")
install.packages('e1071')
library(Biostrings)
library(rentrez)
library(tidyverse)
library(ggplot2)
library(reshape2)
library(data.table)
library(bioseq)
library(randomForest)
library("qpcR")
options(rgl.useNULL = TRUE)
library('rgl')
library('e1071')





#Searches the terms that can be used to extract genes with regards to the database terms. 
entrez_db_searchable("nuccore")


#Extracting HLA-A and B from the database (with web_history)
HLA_A_gene <- entrez_search(db = "nuccore", term = "Human[ORGN] AND  HLA-A[GENE]" ,retmax = 1000,use_history=T)

HLA_B_gene <- entrez_search(db = "nuccore", term = "Human[ORGN] AND  HLA-B[GENE]", retmax = 1000,use_history=T)



####DATA Acquisition for HLA-A sequences

#Getting the sequence of HLA-A with web-history
web_history_HLA_A_seq <- entrez_fetch(db = "nuccore", web_history = HLA_A_gene$web_history, rettype = "fasta")

#Checking the class
class(web_history_HLA_A_seq)


#Saving the HLA-A sequence to the hard drive and separating based on new line 
write(web_history_HLA_A_seq, "HLA_A_seq.fasta", sep = "\n")


# Read it back in as DNA StringSet using the readDNAStringSet() function.
HLA_A_seq_stingset <- readDNAStringSet("HLA_A_seq.fasta")



#Creating a dataframe for the HLA-A
df_HLA_A_seq <- data.frame(HLA_A_Title = names(HLA_A_seq_stingset), HLA_A_Sequence = paste(HLA_A_seq_stingset))
df_HLA_A_seq


#### DATA Exploration for HLA-A


#Calculating the mean length of the sequences for HLA-A, used downstream for creating a barplot
HLA_A_Mean_length = sum(width(HLA_A_seq_stingset)) / length(HLA_A_seq_stingset)
HLA_A_Mean_length

#Summary statistics  of the HLA-A gene sequence 
summary(nchar(HLA_A_seq_stingset))

#checking for NA in HLA-A sequences 
sum(is.na(df_HLA_A_seq$HLA_A_Sequence))


#Getting the range of sequences,this is to know the lower end of sequence length to remove sequences the are short which results in them being compressed when displayed
seq_nchar(dna(df_HLA_A_seq$HLA_A_Sequence)) %>% range()

# Removing sequences that are less than 500bp because some sequences are displayed in a squished manner due to being to small compared to the rest of the sequences
df_HLA_A_seq = df_HLA_A_seq[(which(nchar(df_HLA_A_seq$HLA_A_Sequence)  > 500)),]


#Removing sequences with N character in the sequence 
df_HLA_A_seq <- df_HLA_A_seq[grepl("[^N]*", df_HLA_A_seq$HLA_A_Sequence), , drop = FALSE]

#Calculating the q1 and q3 for HLA-A sequence,this is to filter sequences that are less than and greater than q1 and q3 
HLA_A_Sequence_q1 <- quantile(nchar(df_HLA_A_seq$HLA_A_Sequence), probs = 0.25, na.rm = TRUE)
HLA_A_Sequence_q1

HLA_A_Sequence_q3 <- quantile(nchar(df_HLA_A_seq$HLA_A_Sequence), probs = 0.75, na.rm = TRUE)
HLA_A_Sequence_q3


#Filtering sequences based on q3 and q1 values 
df_HLA_A_seq = df_HLA_A_seq[(which(nchar(df_HLA_A_seq$HLA_A_Sequence)  <= HLA_A_Sequence_q3)),]

df_HLA_A_seq = df_HLA_A_seq[(which(nchar(df_HLA_A_seq$HLA_A_Sequence)  >= HLA_A_Sequence_q1)),]


#Checking that the filtering worked 
summary(nchar(df_HLA_A_seq$HLA_A_Sequence))


####Acquiring HLA-B sequence

#Getting the sequence of HLA-B with web-history
web_history_HLA_B_seq <- entrez_fetch(db = "nuccore", web_history = HLA_B_gene$web_history, rettype = "fasta")

#Checking the class
class(web_history_HLA_B_seq)


#saving the HLA-B sequence to the hard drive and separating based on new line
write(web_history_HLA_B_seq, "HLA_B_seq.fasta", sep = "\n")


# Read it back in as DNA StringSet using the readDNAStringSet() function.
HLA_B_seq_stingset <- readDNAStringSet("HLA_B_seq.fasta")



#Creating a dataframe for the HLA-B
df_HLA_B_seq <- data.frame(HLA_B_Title = names(HLA_B_seq_stingset), HLA_B_Sequence = paste(HLA_B_seq_stingset))
df_HLA_B_seq


#### DATA Exploration for HLA-B

#calculating the mean length of the sequences for HLA-B used downstream to create a barplot 
HLA_B_Mean_length = sum(width(HLA_B_seq_stingset)) / length(HLA_B_seq_stingset)
HLA_B_Mean_length

#Summary statistics  of the HLA-B gene sequence 
summary(nchar(HLA_B_seq_stingset))

#checking for NA in HLA-B sequences 
sum(is.na(df_HLA_B_seq$HLA_B_Sequence))

#Getting the range of sequences,this is to know the lower end of sequence length to remove sequences the are short which results in them being compressed when displayed
seq_nchar(dna(df_HLA_B_seq$HLA_B_Sequence)) %>% range()

# Removing sequences that are less than 500bp because some sequences are displayed in a squished manner due to being to small compared to the rest of the sequences
df_HLA_B_seq = df_HLA_B_seq[(which(nchar(df_HLA_B_seq$HLA_B_Sequence)  > 500)),]


#Removing sequences with N character in the sequence 
df_HLA_B_seq <- df_HLA_B_seq[grepl("[^N]*", df_HLA_B_seq$HLA_B_Sequence), , drop = FALSE]


#Calculating the q1 and q3 for HLA-B sequence,this is to filter sequences that are less than and greater than q1 and q3 respectively 
HLA_B_Sequence_q1 <- quantile(nchar(df_HLA_B_seq$HLA_B_Sequence), probs = 0.25, na.rm = TRUE)
HLA_B_Sequence_q1

HLA_B_Sequence_q3 <- quantile(nchar(df_HLA_B_seq$HLA_B_Sequence), probs = 0.75, na.rm = TRUE)
HLA_B_Sequence_q3


#Filtering sequences based on q3 and q1 values 
df_HLA_B_seq = df_HLA_B_seq[(which(nchar(df_HLA_B_seq$HLA_B_Sequence)  <= HLA_B_Sequence_q3)),]

df_HLA_B_seq = df_HLA_B_seq[(which(nchar(df_HLA_B_seq$HLA_B_Sequence)  >= HLA_B_Sequence_q1)),]

#Checking that the filtering worked 
summary(nchar(df_HLA_B_seq$HLA_B_Sequence))



####Data Analysis 

#calculating the GC frequency in HLA-A and B sequence.Knowing the GC content is important because it causes bias in sequencing data 
df_GC_content <- qpcR:::cbind.na(data.frame(x=letterFrequency(DNAStringSet(df_HLA_B_seq$HLA_B_Sequence), as.prob = TRUE,letters = c("GC"))),data.frame(x=letterFrequency(DNAStringSet(df_HLA_A_seq$HLA_A_Sequence), as.prob = TRUE,letters =c("GC"))))
df_GC_content

#calculating the AT frequency in HLA-A and B sequence.
df_AT_content <- qpcR:::cbind.na(data.frame(x=letterFrequency(DNAStringSet(df_HLA_B_seq$HLA_B_Sequence), as.prob = TRUE,letters = c("AT"))),data.frame(x=letterFrequency(DNAStringSet(df_HLA_A_seq$HLA_A_Sequence), as.prob = TRUE,letters =c("AT"))))
df_AT_content

#column names
Column_names <- c("HLA_B","HLA_A")

#Assigning column names to dataframe
colnames(df_GC_content) <- Column_names
colnames(df_AT_content) <- Column_names

#Adding a Third column to represent the GC content 
df_GC_content <- cbind(Dinucleotide = 'GC', df_GC_content)

df_AT_content <- cbind(Dinucleotide = 'AT', df_AT_content)
df_AT_content



#Merging the GC% dataframe with AT% dataframe to assist in producing the boxplot 
DF_GC_AT_content1 <- rbind(df_GC_content,df_AT_content)


#Changing from long data to wide data to produce the boxplot 
df_long <- reshape2::melt(DF_GC_AT_content1, id.vars=c("Dinucleotide") ,variable.name= "Gene",value.name = "Percentage" ) 

df1_long <- na.omit(df_long)
#Creating the Boxplot
Frequency_Boxplot <- ggplot(df1_long, aes(x=Gene,y=Percentage,fill=Gene))+
  geom_boxplot() + labs(title="AT & GC Frequency") +facet_wrap(~Dinucleotide)+
  theme(plot.title = element_text(hjust = 0.5))

Frequency_Boxplot

#Mean length barplot for HLA-A and B
Mean_len <- data.frame(Gene=c("HLA-A","HLA-B"),Sequence_Mean_length= c(HLA_A_Mean_length,HLA_B_Mean_length) )


#Creating a Mean length bar plot for both genes 
Mean_length_barplot<-ggplot(data=Mean_len, aes(x=Gene, y=Sequence_Mean_length)) +
  geom_bar(stat="identity",fill="steelblue", width=0.5)+
  labs(y= "Sequence Mean Length (bp)", x = "Gene", title="Average Length of HLA-A and HLA-B Sequences")+
  geom_text(aes(label=Sequence_Mean_length), vjust=1.6, color="white",size=3.5)+
  theme(plot.title = element_text(hjust = 0.5))

  
Mean_length_barplot

### Creating the predictors that will be used for RandomForest 

#Creating a dataframe that contains Columns with the Genes and Column with corresponding sequences, This will be used as the training set 
df_HLA_B_seq$Gene <- 'HLA-B'
df_HLA_A_seq$Gene <- 'HLA-A'
df_Combined <- data.table(merge(df_HLA_B_seq,df_HLA_A_seq,by="Gene",all=TRUE))


df_final  <-  melt(df_Combined,
                id.vars="Gene",
                measure.vars=c("HLA_A_Sequence", "HLA_B_Sequence"),
                value.name="Sequence")[order(by=Gene)][, c("Gene", "Sequence")]

#Converting the sequence to DNAstringset class so biostring package is functional on the sequence
df_final <- as.data.frame(df_final)
df_final <- na.omit(df_final)
df_final$Sequence <- DNAStringSet(df_final$Sequence)

#Calculating the nucleotide frequencies and appending onto our dataframe using cbind()
df_final <- cbind(df_final, as.data.frame(letterFrequency(df_final$Sequence, letters = c("A", "C","G", "T"))))


#Adding A, T, and G proportions in relation to total nucleotides 
df_final$Aprop <- (df_final$A) / (df_final$A + df_final$T + df_final$C + df_final$G)

df_final$Tprop <- (df_final$T) / (df_final$A + df_final$T + df_final$C + df_final$G)

df_final$Gprop <- (df_final$G) / (df_final$A + df_final$T + df_final$C + df_final$G)


#Adding dinucleotide frequency (k-mers of length 2), here using proportions as that let's us account for variability in sequence lengths.
df_final <- cbind(df_final, as.data.frame(dinucleotideFrequency(df_final$Sequence, as.prob = TRUE)))


#Adding trinucleotide frequency (k-mers of length 3)
df_final <- cbind(df_final, as.data.frame(trinucleotideFrequency(df_final$Sequence, as.prob = TRUE)))
ncol(df_final)

#####TRAINING CLASSIFICATION MODEL

#converting string format back to character data, so we can use tidyverse functions.
df_final$Sequence <- as.character(df_final$Sequence)

#Creating a Validation set and a Training set which don't overlap, The Training set will be used to train the algorithm and the validation set will be used for the prediction. Set.seed is used for reproducible.
set.seed(200)

dfValidation <- df_final %>%
  group_by(Gene) %>%
  sample_n(1000)



set.seed(13)
df_Training <- df_final %>%
  filter(!Sequence %in% dfValidation$Sequence) %>%
  group_by(Gene) %>%
  sample_n(61)

###RandomForest
#Training the classifier using the training set Using the randomforest Algorithm 
gene_classifier <- randomForest::randomForest(x = df_Training[, 3:89], y = as.factor(df_Training$Gene), ntree = 50, importance = TRUE)

gene_classifier

#Using the classifier and the validation set to separate HLA-A and HLA-B gene
predictValidation <- predict(gene_classifier, dfValidation[, 3:89])

predictValidation

#Create our own confusion matrix.
confusion_matrix <- table(observed = dfValidation$Gene, predicted = predictValidation)
confusion_matrix


#Creating a data frame with Randomforest confusion matrix.
CM_RF_table <- data.frame(confusion_matrix)

Column_names_2 <- c("Observed","Predicted","Freq")
colnames(CM_RF_table) <- Column_names_2

#Creating a plot for the RF data 
RF_plot<- ggplot(data =  CM_RF_table, mapping = aes(x = Predicted, y = Observed)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none")+
  labs(title = "Random Forest")+
  theme(plot.title = element_text(hjust = 0.5))

RF_plot



###Naive Bayes 
#Training the classifier using the training set using naive Bayes algorithm 

gene_classifier_Bayes <- naiveBayes(df_Training[, 3:89], df_Training$Gene)

gene_classifier_Bayes

#Using the classifier and the validation set to separate HLA-A and HLA-B gene

predictValidation <- predict(gene_classifier_Bayes, newdata=dfValidation[, 3:89])

predictValidation

#Create our own confusion matrix.
cm_Bayes<- table(dfValidation$Gene, predictValidation)

cm_Bayes


#Creating a data frame with Naive bayes confusion matrix.
cm_Bayes_table <- data.frame(cm_Bayes)

Column_names_2 <- c("Observed","Predicted","Freq")
colnames(cm_Bayes_table) <- Column_names_2

#Creating a plot for the Naive bayes data 
Bayes_plot<- ggplot(data =  cm_Bayes_table, mapping = aes(x = Predicted, y = Observed)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() + theme(legend.position = "none")+
  labs(title = "Naive Bayes")+
  theme(plot.title = element_text(hjust = 0.5))

Bayes_plot


knitr::stitch('Assignment2.r')
