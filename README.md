[Uploading cavity_flow_readme.md…]()
# Lid-Driven Cavity Flow Simulation

A high-performance OpenMP-parallelized Computational Fluid Dynamics (CFD) solver for the classic lid-driven cavity flow problem using the vorticity-streamfunction formulation.

## Table of Contents
- [Problem Description](#problem-description)
- [Mathematical Formulation](#mathematical-formulation)
- [Numerical Method](#numerical-method)
- [Features](#features)
- [Requirements](#requirements)
- [Compilation and Execution](#compilation-and-execution)
- [Input Parameters](#input-parameters)
- [Output Files](#output-files)
- [Performance Considerations](#performance-considerations)
- [Validation](#validation)
- [Troubleshooting](#troubleshooting)
- [References](#references)

## Problem Description

The lid-driven cavity flow is a fundamental benchmark problem in computational fluid dynamics. It consists of a square cavity filled with viscous fluid, where the top wall moves with a constant velocity while all other walls remain stationary. This creates a recirculating flow pattern with primary and secondary vortices.

### Physical Setup
- **Domain**: Unit square cavity [0,1] × [0,1]
- **Boundary Conditions**:
  - Top wall (lid): Moving with velocity u = 1, v = 0
  - Bottom, left, right walls: No-slip conditions (u = v = 0)
- **Fluid Properties**: Incompressible, Newtonian fluid

## Mathematical Formulation

### Governing Equations
The incompressible Navier-Stokes equations are solved using the vorticity-streamfunction (ω-ψ) formulation:

#### Vorticity Transport Equation
```
∂ω/∂t + u·∇ω = (1/Re)∇²ω
```

#### Streamfunction Poisson Equation
```
∇²ψ = -ω
```

#### Velocity-Streamfunction Relations
```
u = ∂ψ/∂y,  v = -∂ψ/∂x
```

where:
- ω = vorticity (scalar in 2D)
- ψ = streamfunction
- Re = Reynolds number = UL/ν
- U = characteristic velocity (lid velocity = 1)
- L = characteristic length (cavity width = 1)
- ν = kinematic viscosity

### Boundary Conditions

#### Streamfunction (Dirichlet)
- ψ = 0 on all walls (constant streamfunction along solid boundaries)

#### Vorticity (Computed from no-slip condition)
- **Stationary walls**: ω = 2(ψ_wall - ψ_interior)/Δn²
- **Moving top lid**: ω = 2(ψ_wall - ψ_interior)/Δy² - 2U_lid/Δy

## Numerical Method

### Discretization
- **Spatial**: Second-order central finite differences on uniform Cartesian grid
- **Temporal**: Steady-state solution using iterative methods
- **Grid**: 128×128 uniform mesh (customizable)

### Solution Algorithm
1. **Vorticity Update**: Solve transport equation using upwind-differenced convection terms
2. **Streamfunction Update**: Solve Poisson equation using 5-point stencil
3. **Boundary Update**: Compute wall vorticity from current streamfunction
4. **Convergence Check**: Monitor L2 norm of field changes

### Iterative Method
- **Jacobi Iteration**: Thread-safe parallel implementation
- **Under-relaxation**: Improves convergence stability (α = 0.7)
- **Convergence Criterion**: RMS error < 1×10⁻⁶ for both ω and ψ

### Numerical Stability
- **Upwind Convection**: Prevents numerical instabilities at high Reynolds numbers
- **Under-relaxation**: Ensures convergence for challenging flow conditions
- **Grid Resolution**: Adequate mesh density to resolve boundary layers

## Features

### High Performance Computing
- **OpenMP Parallelization**: Efficient multi-threading for all major computational loops
- **Thread-Safe Design**: Jacobi iterations with separate read/write buffers
- **NUMA Awareness**: Optimized memory access patterns
- **Scalable Architecture**: Performance scales with available CPU cores

### Robustness
- **Stable Numerics**: Upwind differencing and under-relaxation for stability
- **Convergence Monitoring**: Real-time progress tracking and error reporting
- **Automatic Termination**: Maximum iteration limits prevent infinite loops
- **Error Handling**: Robust file I/O and memory management

### Comprehensive Output
- **Multiple Data Formats**: Tecplot-compatible output files
- **Field Visualization**: Complete streamfunction and vorticity fields
- **Profile Data**: Centerline velocity profiles for validation
- **Progress Tracking**: Iteration-by-iteration convergence monitoring

## Requirements

### System Requirements
- **Operating System**: Linux, macOS, or Windows with appropriate compiler
- **Memory**: Minimum 1 GB RAM (for default 128×128 grid)
- **CPU**: Multi-core processor recommended for optimal performance

### Software Dependencies
- **C Compiler**: GCC 4.9+, Intel ICC, or Clang with OpenMP support
- **OpenMP**: Version 3.0 or higher
- **Math Library**: Standard math library (libm)

### Optional Tools
- **Visualization**: Tecplot, ParaView, or matplotlib for results visualization
- **Profiling**: Intel VTune, GNU gprof, or similar performance analysis tools

## Compilation and Execution

### Basic Compilation
```bash
gcc -O3 -fopenmp -lm cavity_flow.c -o cavity_flow
```

### Optimized Compilation (Intel Compiler)
```bash
icc -O3 -qopenmp -xHost cavity_flow.c -o cavity_flow
```

### Compilation Flags Explained
- `-O3`: Maximum optimization level
- `-fopenmp` / `-qopenmp`: Enable OpenMP support
- `-lm`: Link math library
- `-xHost` (Intel): Optimize for host processor architecture

### Execution
```bash
# Set number of OpenMP threads
export OMP_NUM_THREADS=4

# Run the simulation
./cavity_flow
```

### Thread Configuration
```bash
# Use all available cores
export OMP_NUM_THREADS=$(nproc)

# Set thread affinity for better performance
export OMP_PROC_BIND=true
export OMP_PLACES=cores
```

## Input Parameters

### Primary Parameters (modify in source code)
```c
m = 128;              // Grid points in x-direction
n = 128;              // Grid points in y-direction
Re = 400.0;           // Reynolds number
alpha_psi = 0.7;      // Streamfunction under-relaxation
alpha_omega = 0.7;    // Vorticity under-relaxation
max_iter = 20000;     // Maximum iterations
```

### Parameter Guidelines

#### Grid Resolution
- **Coarse (64×64)**: Fast computation, lower accuracy
- **Standard (128×128)**: Good balance of speed and accuracy
- **Fine (256×256)**: High accuracy, longer computation time
- **Very Fine (512×512)**: Research-grade accuracy, substantial computational cost

#### Reynolds Number Effects
- **Re < 100**: Steady, simple flow pattern
- **Re = 100-1000**: Development of secondary vortices
- **Re = 1000-5000**: Complex vortex structures
- **Re > 5000**: Requires very fine grids and special treatment

#### Under-relaxation Values
- **α = 0.5**: Very stable, slow convergence
- **α = 0.7**: Good balance (recommended)
- **α = 0.9**: Fast convergence, may be unstable
- **α = 1.0**: No relaxation, often unstable

## Output Files

### Generated Files
1. **Stream_cavity_flow_Re400.dat**: Complete streamfunction field
2. **Vorticity_cavity_flow_Re400.dat**: Complete vorticity field
3. **u_at_centreline_Re400.dat**: u-velocity along vertical centerline
4. **v_at_centreline_Re400.dat**: v-velocity along horizontal centerline

### File Format
All files use Tecplot ASCII format with headers:
```
VARIABLES = X, Y, PSI
ZONE I=128, J=128, F=POINT
x_1  y_1  value_1
x_2  y_2  value_2
...
```

### Data Interpretation
- **Streamfunction**: Contour lines represent flow streamlines
- **Vorticity**: Positive values indicate counterclockwise rotation
- **Velocity Profiles**: Used for quantitative validation against literature

## Performance Considerations

### Optimization Strategies

#### Thread Scaling
```bash
# Benchmark different thread counts
for threads in 1 2 4 8 16; do
    export OMP_NUM_THREADS=$threads
    time ./cavity_flow
done
```

#### Memory Optimization
- Use static arrays to avoid stack overflow
- Consider memory bandwidth limitations
- NUMA-aware memory allocation for large systems

#### Computational Optimization
- Compiler vectorization with `-O3` flag
- Processor-specific optimizations (`-march=native`)
- Profile-guided optimization (PGO) for production runs

### Performance Expectations

#### Typical Performance (128×128 grid)
- **Single Thread**: ~30-60 seconds
- **4 Threads**: ~10-20 seconds
- **8 Threads**: ~5-15 seconds
- **16+ Threads**: Diminishing returns due to memory bandwidth

#### Scaling Characteristics
- **Good Scaling**: Up to number of physical cores
- **Memory Bound**: Performance limited by bandwidth at high core counts
- **NUMA Effects**: Performance drops on multi-socket systems without optimization

## Validation

### Reference Solutions
Compare results with established benchmarks from literature:
- Ghia et al. (1982) - Standard benchmark data
- Bruneau & Saad (2006) - High Reynolds number results
- Erturk et al. (2005) - Fine-grid reference solutions

### Validation Metrics
1. **Centerline Velocity Profiles**: Compare u(0.5,y) and v(x,0.5)
2. **Vorticity Distribution**: Check primary vortex location and strength
3. **Streamfunction Contours**: Verify flow pattern characteristics
4. **Global Flow Properties**: Mass conservation and energy balance

### Expected Results (Re = 400)
- **Primary Vortex Center**: ~(0.56, 0.61)
- **Maximum Streamfunction**: ~0.114
- **Secondary Vortices**: Small corner recirculation zones

## Troubleshooting

### Common Issues

#### Compilation Problems
```bash
# Missing OpenMP support
sudo apt-get install libomp-dev  # Ubuntu/Debian
brew install llvm                # macOS with Homebrew

# Missing math library
gcc -O3 -fopenmp cavity_flow.c -lm -o cavity_flow
```

#### Runtime Issues
```bash
# Segmentation fault - increase stack size
ulimit -s unlimited

# Slow convergence - adjust under-relaxation
# Modify alpha_psi and alpha_omega in source code
```

#### Performance Problems
```bash
# Check thread utilization
htop  # Monitor CPU usage during execution

# Verify OpenMP is working
export OMP_DISPLAY_ENV=true
./cavity_flow
```

### Convergence Issues
- **Oscillating residuals**: Reduce under-relaxation parameters
- **Slow convergence**: Check grid quality and boundary conditions
- **Divergence**: Lower Reynolds number or refine grid

### Memory Issues
- **Stack overflow**: Use `ulimit -s unlimited` or make arrays static
- **Insufficient RAM**: Reduce grid size or use out-of-core algorithms

## References

### Key Literature
1. **Ghia, U., Ghia, K. N., & Shin, C. T. (1982)**  
   "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method"  
   *Journal of Computational Physics*, 48(3), 387-411

2. **Bruneau, C. H., & Saad, M. (2006)**  
   "The 2D lid-driven cavity problem revisited"  
   *Computers & Fluids*, 35(3), 326-348

3.
