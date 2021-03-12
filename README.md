# The Piece-Wise Spline Wigner-Ville Distribution MATLAB Package V1.0
This repository supplements the following publication:
M. Al-Sa'd, B. Boashash, and M. Gabbouj, “Design of an Optimal Piece-Wise Spline Wigner-Ville Distribution for TFD Performance Evaluation and Comparison”, IEEE Transactions on Signal Processing, 2021

It contains MATLAB scripts, functions, and datasets to produce results, and supporting figures illustrated in the paper.
Demo Scripts:
The developed PW-WVD package contains the following demo scripts within its directory:
1.	Demo_1_example_simulation.m
    •	Description: This demo script produces the results that are depicted in Fig. 2 of the paper.
    •	Process: It generates an example 2-component non-stationary signal and computes its PW-WVD, WVD, and the WVD cross-terms. In addition, it shows the signal IA, IF, and time-support functions.
2.	Demo_2_tfd_comparison.m
•	Description: This demo script produces the results that are depicted in Fig. 3 of the paper.
•	Process: It generates an example 2-component non-stationary signal and computes the signal ideal TFD, PW-WVD, and the MPD using chirplet atoms. Note that the MPD is precomputed and saved in MPD_example.mat.
3.	Demo_3_resolution_visualization.m
•	Description: This demo script produces the illustration in Fig. 4 of the paper.
•	Process: It yields a one-dimensional visual interpretation for the proposed TFD resolution measure.
