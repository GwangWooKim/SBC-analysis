<img src="/imgs/SBC.jpg" width="50%" height="50%">
Myungin Baek et al., Molecular Logic of Spinocerebellar Tract Neuron Diversity and Connectivity, Cell Reports, 2019.

# SBC-analysis
This repository contains codes which are used in analyzing the data from DISIST. Note that 
* We can't provide the used data (`merged_count_matrix_3`) and its metadata (`condtion_3`)
* During analyzing, __we figure out that the used dataset has problems (in preprocessing-level). So DIGIST has been resolving this problem. We will re-analyze the corrected data and report results.__
* The code 92-113 might not be able to work. I will deal with this part after getting the data.

## Background
Suppose we have a `n_gene` by `n_cell` matrix. (it is ok to understand a gene as a feature and a cell as one datapoint if you want). Each cell follows from a known population. In this situation, which genes represent a populations? how much represent? 

## Method
Statistical models give us good answers. "How much" is arranged by p-value and "Which" is determined by genes having lower p-value. But, if there is [covariate effect](https://en.wikipedia.org/wiki/Batch_effect), we should adjust the original data or establish a corrected data matrix. Our data contains two different sequencing data. One is from DropSeq and the other is from Smart-seq2. They are much different. One way to adress this issue is to introduce a surrogate variable, which will be used at modeling step with the data. Please refer to [`Draft_1_4.pptx`](https://github.com/GwangWooKim/SBC-analysis/blob/main/Draft_1_4.pptx) in this repository if you want to know more backgorunds, methods, and results. (NB: the defintion of our surrogate variable is still private.)

## Results (partial)
<img src="/imgs/result1.png" width="70%" height="70%">
This figure shows that 40% genes are filtered out (that is, they are not tested whether each of them is meaningful or not)
<img src="/imgs/result2.png" width="70%" height="70%">
The same result of the previous figure, but different aspect. Indeed, the more meaningless (higher p-value) it was, the less tested it was.
<img src="/imgs/result3.png" width="100%" height="100%">
A final output of DESeq2 with some genes
<img src="/imgs/result4.png" width="100%" height="100%">
A gene clustering result. More results can be fonund in Draft_1_4.pptx.
