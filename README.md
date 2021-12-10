# LOA_correlation
This is the code that was used to obtain the curves for Fig.4a in the following paper: "Effect of degree correlations above the first shell on the percolation transition" Valdez et al. EPL (Europhysics Letters) 96 (3), 38001 (2011)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Input:

N_node: number of nodes

election: 1 (if the user wants an assortative network) or -1 (if the user wants a disassortative network)

N_iterat: number of iterations to correlate the network. The Pearson correlation coefficient increases with "N_iterat"

nrea: number of realizations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Output:

Pinf*dat: fraction of nodes belonging to the giant compoenent as a function of p for a link percolation process

Smed*.dat: mean finite cluster size as a function of p

PearsonCorrCoef*.dat: average Pearson correlation coefficient

SmedMax*.dat": mean height of the peak of for a network with "N_node" nodes.
