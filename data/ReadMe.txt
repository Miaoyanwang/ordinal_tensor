###########################
InCar_Music.RData

InCarMusic is a mobile application (Android) that offers music recommendations to passengers based on contexts.

The R object ``tensor'' records the ratings (on 1-5 scale) of 42 users to 139 songs on 26 contexts. The ``context'' attributes describe the 26 contextual conditions. The ``music'' attributes describe the 139 songs and their corresponding genres.

Two goals:
1. tensor completion: offer context-specific music recommendations to users. 
2. dimension reduction: identify similar songs and similar contexts based on rating pattern. 

Reference: Baltrunas, Linas, et al. "Incarmusic: Context-aware music recommendations in a car." E-Commerce and Web Technologies. Springer Berlin Heidelberg, 2011. 89-100.

###########################
dti_brain.RData
The dti tensor is a 68 × 68 × 136 ordinal tensor consisting of structural connectivity patterns among 68 brain regions for 136 individuals from Human Connectome Project (HCP). All the individual images were preprocessed following a standard pipeline (Zhang et al., 2018), and the brain was parcellated to 68 regions-of-interest following the Desikan atlas (Desikan et al., 2006). The tensor entries encode the strength (on 1-3 scale) of fiber connections between 68 brain regions for each of the 136 individuals.

Goal: 
Dimension reduction: identify similar brain nodes based on connectivity pattern. 

Reference: Section 7.1 https://arxiv.org/pdf/1910.09499.pdf

###########################
BrainNet is a useful software for visualizing brain connections. 

Instructions (The detailed instruction is brainNet_Manual.pdf):
1. Install Matlab (free) from UW-Madison
2. Open matlab and type "BrainNet"
3. Load files:
Surface file: BrainNetView.../Data/miaoyan/BrainMesh_ICBM152.nv

Data file (node): BrainNetView.../Data/miaoyan/Desikan_miaoyan.node.txt 

Data file (edge): A user-specified input. This is a matrix file representing the brain connection. The edge can be raw connection strength, estimated strength, truncated strength, or any measures you want to plot over the brain edges. I included example0.txt file as an example input.

Note: (You may need enable "all files *" option by clicking the "options" in the pop-up console.)

4. Click okay and play around with various settings.