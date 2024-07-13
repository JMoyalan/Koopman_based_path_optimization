<p align="center">
<img src="figures/off_road_navigation.PNG" width="1000">
</p>

## Overview
The code in this repository considers the problem of optimal control design for navigation on off-road terrain. We use a traversability measure to characterize the difficulty of navigation on off-road terrain. The traversability measure captures terrain properties essential for navigation, such as elevation maps, roughness, slope, and texture. We provide a convex formulation to the off-road navigation problem by lifting the problem to the density space using the linear Koopman and Perron-Frobenius (P-F) operators. Our convex framework is then implemented for path planning of the legged robots in an unstructured environment. Please look at the following papers for high-level details of the framework.
```
1. Moyalan, Joseph, Yongxin Chen, and Umesh Vaidya. "Convex approach to data-driven off-road navigation via linear transfer operators." IEEE Robotics and Automation Letters 8.6 (2023): 3278-3285.
2. Moyalan, Joseph, et al. "Off-Road Navigation of Legged Robots Using Linear Transfer Operators." IFAC-PapersOnLine 56.3 (2023): 613-618.
```
