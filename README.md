# GKLOMLI

Paper: GKLOMLI: A link prediction model for inferring miRNA-lncRNA interactions by using Gaussian kernel-based method on network profile and linear optimization algorithm
Author: Leon Wong, Lei Wang, Zhu-Hong You, Chang-An Yuan, Yu-An Huang, and Mei-Yuan Cao

****
# File description
****

## Dataset:
### Expression_profile_data: 
#### invalid_association_expression.mat: The ids of known miRNA-lncRNA interactions without either of lncRNA and miRNA expression profiles.
#### invalid_lnc_expression.mat: The ids of lncRNAs without lncRNA expression profiles.
#### invalid_mi_expression.mat: The ids of lncRNAs without lncRNA expression profiles.
#### lnc_expression_similarity_matrix.mat: The lncRNA-lncRNA similarity matrix of expression profiles.
#### mi_expression_similarity_matrix.mat: The miRNA-miRNA similarity matrix of expression profiles.
----------------------------
## Function_data:
### invalid_association_function.mat: The ids of known miRNA-lncRNA interactions without either of lncRNA and miRNA function feature.
### invalid_lnc_function.mat: The ids of lncRNAs without lncRNA function feature.
### invalid_mi_function.mat: The ids of lncRNAs without lncRNA function feature.
### lnc_function_similarity_matrix.mat: The lncRNA-lncRNA similarity matrix of function.
### mi_function_similarity_matrix.mat: The miRNA-miRNA similarity matrix of function.

## Sequence_data:
### invalid_lnc_seq.mat: The ids of lncRNAs without lncRNA sequence information.
### lnc_seq_similarity_matrix.mat: The lncRNA-lncRNA similarity matrix of sequence.
### mi_seq_similarity_matrix.mat: The miRNA-miRNA similarity matrix of sequence.
----------------------------
## miRNA-lncRNA_interaction:
### lncRNA_id.txt The id list of lncRNA name.
### miRNA_id.txt The id list of miRNA name.
### miRNA-lncRNA_idlist.txt: The miRNA-lncRNA interaction data downloaded from lncRNASNP.
### interaction.mat: The miRNA-lncRNA interaction network is represented as an adjacent matrix.
****
## Tools:
### gaussiansimilarity.m: The method is to calculate the similarity based on Gaussian kernel in term of interaction profile.
### Normalize.m: The method is to normalize all elements with a range from 0 to 1.  
### positiontooverallauc.m: The codes for plotting the ROC curves based on the ranks of testing samples.
### Method_Gaussian.m:  Link prediction method by using Gaussian kernel-based method on network profile for the implemention of k-fold cross validation.
### Method_SP.m: Link prediction method by using similarity measurement methods on different type of bio-profiles for the implemention of k-fold cross validation.
### Method_Single.m: Link prediction method without any side information for the implemention of k-fold cross validation.
****
## main_expression.m: Tha main function for inferring miRNA-lncRNA interactions by using expression profile.
## main_function.m: Tha main function for inferring miRNA-lncRNA interactions by using functional profile.
## main_seq.m: Tha main function for inferring miRNA-lncRNA interactions by using sequence profile.
## main_noProfile.m: Tha main function for inferring miRNA-lncRNA interactions without side information.
## main_expression.m: Tha main function for inferring miRNA-lncRNA interactions by using network similarity.
