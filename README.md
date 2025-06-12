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
- **`GBM.m`**  
  Compares **weak convergence** of Euler-Maruyama discretization against exact solutions and against a finer time-step (2h method). Demonstrates error reduction with finer timesteps.

### 4. Survival Probability Estimation in Population Models
- **`malthusian.m`**  
  Estimates survival probability in a Malthusian growth model using **log-normal GBM**. Incorporates up-in and down-out barrier conditions. Calculates probability, error, and computation time.

---

## 🔧 Requirements

- MATLAB R2020+ (or compatible)
- No external toolboxes required
- Sufficient memory for large simulations (1e7 samples used in some files)

---

## ▶️ How to Run

Each script is a standalone function. You can call them with appropriate parameters. For example:

```matlab
[V, ster, CPUt, varsc, eb] = ouscratch_bb(1e5, 64, 42, 1);
[prob, err, ~, time, var, paths] = malthusian(1e6, 100, 1234);
