# Stochastic Simulation and Barrier Option Pricing (MATLAB)

This repository contains a set of MATLAB scripts for simulating stochastic differential equations (SDEs), evaluating barrier options, and estimating survival probabilities in population models using Monte Carlo methods. The focus is on comparing different techniques for dealing with discrete-time approximations and barrier-crossing events (e.g., Brownian Bridge, barrier shifting).

---

## 📁 Contents

### 1. Option Pricing (Ornstein–Uhlenbeck Process)
- **`ouscratch_bb.m`**  
  Simulates barrier options under an Ornstein–Uhlenbeck (OU) process using **Brownian Bridge correction** to detect barrier crossings within time steps.

- **`ouscratch_shifted.m`**  
  Simulates the same setup using the **barrier shifting method** as a bias correction instead of Brownian Bridge.

### 2. Brownian Bridge Utilities
- **`brownianBridge.m`**  
  Generates sample paths of a Brownian Bridge. Used by the OU pricing models to interpolate barrier crossings.

### 3. Geometric Brownian Motion (GBM) Weak Convergence Study
- **`GBM_WeakConvergence_Simulation.m`**  
  Compares **weak convergence** of Euler-Maruyama discretization against exact solutions and against a finer time-step (2h method). Demonstrates error reduction with finer timesteps.

### 4. Survival Probability Estimation in Population Models
- **`malthusian.m`**  
  Estimates survival probability in a Malthusian growth model using **log-normal GBM**. Incorporates up-in and down-out barrier conditions. Calculates probability, error, and computation time.

---

### 📄 Report Summary

The accompanying PDF file `BARrier_options.pdf` contains the full report for this project. The primary objective of this practical work is to investigate the numerical challenges in simulating **path-dependent options**, particularly focusing on **barrier options**. The report addresses several key goals:

1. Deriving and solving boundary value problems for down-and-out European call options, comparing these with vanilla calls.
2. Implementing **Euler–Maruyama simulations** to study weak convergence in both Gaussian and non-Gaussian settings.
3. Applying techniques like the **Brownian bridge correction** and **barrier shifting** to improve simulation accuracy.
4. Exploring stochastic models such as the **Ornstein–Uhlenbeck process** and **Geometric Brownian Motion (GBM)** in the presence of barriers.
5. Estimating **survival probabilities** in stochastic population models using numerical and semi-analytical approaches.

Through these approaches, the report evaluates bias, statistical error, and computational efficiency to identify effective strategies for accurate and efficient simulation in financial mathematics.


