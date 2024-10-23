.
├── Analysis_Enterotype.xlsx ## A cross-table of which Taxa are significantly different for Bact_bact vs the ones significant in individual that are Bact2 and transition to "other" Enterotype

├── TOP_10_most_different_Bact2_other.csv ## Top 10 edges with the biggest difference-in-difference between Bact2-Bact2 and Bact2-others

├── graph_bact2_bact2.csv # For each edge the value in Bact2-Bact2 individual before and after treatment (w0 vs w24)
├── graph_bact2_other.csv # Same for Bact2-Other individuals
├── graph_other_other.csv # Same for Other-Other individuals

├── graph_diff_bact2_bact2.csv # Difference of the edge values in Bact2-Bact2 individual
├── graph_diff_bact2_other.csv # Same for Bact2-Other individuals
├── graph_diff_other_other.csv # Same for Other-Other individuals

├── node_dist_bact2_bact2.csv # P-value and ranking sum of the distance between all edges belonging to a node. The lower the p-value the most different the edges of that node are before/after treatment in the Bact2-Bact2 individuals 
├── node_dist_bact2_other.csv # Same for Bact2-Other individuals
├── node_dist_other_other.csv # Same for Other-Other individuals

├── significant_nodes_bact2_bact2.txt ## only significant nodes extracted. 
├── significant_nodes_bact2_other.txt
├── significant_nodes_other_other.txt

0 directories, 15 files