<p align="center">
<img src="figures/off_road_navigation.PNG" width="1000">
</p>

## Overview
The code in this repository considers the problem of optimal control design for navigation on off-road terrain. We use a traversability measure to characterize the difficulty of navigation on off-road terrain. The traversability measure captures terrain properties essential for navigation, such as elevation maps, roughness, slope, and texture. We provide a convex formulation to the off-road navigation problem by lifting the problem to the density space using the linear Koopman and Perron-Frobenius (P-F) operators. The implementation of our convex framework is done for path planning of the legged robots in unstructured environment.
