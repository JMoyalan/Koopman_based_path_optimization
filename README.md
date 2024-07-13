<p align="center">
<img src="figures/off_road_navigation.PNG" width="1000">
</p>

## Overview
The code in this repository considers the problem of optimal control design for navigation on off-road terrain. We use a traversability measure to characterize the difficulty of navigation on off-road terrain. The traversability measure captures terrain properties essential for navigation, such as elevation maps, roughness, slope, and texture. We provide a convex formulation to the off-road navigation problem by lifting the problem to the density space using the linear Koopman and Perron-Frobenius (P-F) operators. Our convex framework is then implemented for path planning of the legged robots in an unstructured environment. Please look at [paper 1](https://ieeexplore.ieee.org/abstract/document/10081458) and [paper 2](https://www.sciencedirect.com/science/article/pii/S2405896323024254) for high-level details of the framework.

