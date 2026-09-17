# vortRing

Simulation of axisymmetric, buoyant vortex ring formation from an initial, uniform buoyancy anomaly in the $\mathrm{Re} \to \infty$ and $\mathrm{Pe} \to \infty$ limit.

- Boussinesq approximation for buoyancy
- interfacial vortex sheet between uniform buoyant region and neutral ambient
- vortex sheet represented as a graph to handle topology optimization by surgery
- integrated with 4th-order Runge-Kutta

## regularization
- damping shortest wavelengths in the induced velocity field (Krasny 1986)

## tree algorithm for velocity evaluation
- extending Dynnikova's (2009) method to handle axisymmetric cases
- reduces complexity from $\mathcal{O}(N^2)$ to $\mathcal{O}(N \log N)$

## dynamic rediscretization
- inserting new nodes in stretched regions ($\delta$-adjusted)

## surgery
- reduction of lamellar structures by merging close and parallel segments (topology changes)
- conservative in terms of total circulation $\Gamma$ and $d\Gamma/dt$
- reduces exponential node growth to approximately $\mathcal{O}(N^2)$
- optimized using tree clustering to $\mathcal{O}(N \log N)$
- inspired by Dritschel's (1988) method for isocontours

## parallelization
- hybrid OpenMP + MPI parallelization in a simple fork model with centralized memory on the root and broadcast synchronization

## applications
(Jędrejko 2025)
- hierarchical Kelvin–Helmholtz instability with wavelength doubling
- estimation of the interface fractal dimension

## Build & Execution

### Prerequisites
- C++17 compliant compiler (`g++` or `clang++`)
- CMake 3.14+
- OpenMP
- MPI implementation (e.g., OpenMPI or MPICH)

### Building
```bash
mkdir -p build results
cd build
cmake ..
make -j
```

## references

* **Krasny, R.** (1986). *Desingularization of periodic vortex sheet roll-up*. Journal of Computational Physics, 65(2), 292–313.  
  [https://doi.org/10.1016/0021-9991(86)90210-X](https://doi.org/10.1016/0021-9991(86)90210-X)

* **Dritschel, D. G.** (1988). *Contour surgery: A topological reconnection scheme for extended integrations using contour dynamics*. Journal of Computational Physics, 77(1), 240–266.  
  [https://doi.org/10.1016/0021-9991(88)90165-9](https://doi.org/10.1016/0021-9991(88)90165-9)

* **Dynnikova, G. Ya.** (2009). *Fast technique for solving the $N$-body problem in flow simulation by vortex methods*. Computational Mathematics and Mathematical Physics, 49(8), 1389–1396.  
  [https://doi.org/10.1134/S0965542509080090](https://doi.org/10.1134/S0965542509080090)

* **Jędrejko, P.** (2025). *Turbulent coherent structures in thermal vortex rings*. Journal of Theoretical and Applied Mechanics, 63(4), 903–913.  
  [https://doi.org/10.15632/jtam-pl/208833](https://doi.org/10.15632/jtam-pl/208833)