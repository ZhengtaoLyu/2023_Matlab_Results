# 🧠 Biological Population Dynamics Simulation and Predator-Prey Model

This repository contains **three MATLAB simulations** for modeling:

1. **Stage-structured population growth** (matrix model)
2. **2D heat diffusion (explicit finite difference method)**
3. **Predator-prey swarm dynamics (multi-agent simulation)**

Each part demonstrates different computational modeling approaches in ecology and biological systems.

---

## 📁 Structure Overview

| Code File | Description | Key Techniques |
|-----------|-------------|----------------|
| `code1_population_matrix.m` | Matrix-based projection of insect population stages over time | Linear algebra, dominant eigenvalue analysis |
| `code2_heat_equation.m`     | Solves 2D heat equation using explicit finite difference | Partial differential equations (PDEs), grid discretization |
| `code3_predator_prey_simulation.m` | Multi-agent simulation of predator-prey dynamics with social forces | Swarm behavior, vector math, emergent systems |

---

## 👤 Target Audience

This project is designed for:

- **Technical hiring managers** evaluating analytical modeling and simulation proficiency
- **Engineering leads / data science supervisors** reviewing mathematical modeling, agent-based simulation, and code structure

---

## 🔬 Methodological Background

| Model Type | Algorithm / Method | Reference |
|------------|--------------------|-----------|
| Stage-structured matrix population model | Leslie matrix / transition matrix model | Caswell, H. (2001). *Matrix Population Models: Construction, Analysis, and Interpretation*. |
| Explicit finite difference for 2D heat equation | Central difference in space + forward Euler in time | Smith, G. D. (1985). *Numerical Solution of Partial Differential Equations: Finite Difference Methods*. |
| Predator-prey multi-agent simulation | Boid-style interaction with repulsion/orientation/attraction | Couzin, I. D. et al. (2002). *Collective memory and spatial sorting in animal groups*. *Journal of Theoretical Biology* |

---

## 📊 Technical Highlights

### 1. Stage-structured population growth (`code1_population_matrix.m`)

- Uses a 6×6 transition matrix to simulate life stages: Eggs → Larvae → Pupae → Adults.
- Eigenvalue λ > 1 implies exponential growth.
- Custom function `tfa4management()` tests matrix perturbations for management (e.g., pest suppression).
- Visualization: time-series plots of population dynamics per stage.

📌 *Improvement opportunities*:  
– Automate matrix tuning to hit target eigenvalue  
– Add uncertainty or stochasticity for realistic variability  
– Allow matrix input via file (CSV/JSON) for flexibility

---

### 2. 2D Heat Diffusion (`code2_heat_equation.m`)

- Explicit finite difference scheme for heat diffusion over a 2D grid.
- Boundary initialization with `sin` and `cos` perturbations.
- Stable under CFL condition: spatial and time steps well-controlled.

📌 *Improvement opportunities*:  
– Replace explicit scheme with implicit (e.g., Crank–Nicolson) for better numerical stability  
– Add real-time visualization during computation  
– Vectorize inner loops for speedup

---

### 3. Predator-Prey Simulation (`code3_predator_prey_simulation.m`)

- 2D swarm simulation with 20 prey and 3 predators.
- Prey behavior: repulsion (Rr), orientation (Ro), attraction (Ra).
- Predator behavior: seeks nearest visible prey, with directional noise.
- Emergent flocking and pursuit behavior.
- Animation rendered with `quiver` arrows (velocity vectors).

📌 *Improvement opportunities*:  
– Modularize prey and predator logic into functions/classes  
– Add real-time metrics (prey count, energy loss, clustering index)  
– Introduce boundary conditions or habitat zones

---

## 📽️ Demo Snapshots

### Predator-Prey Flocking System
![Animation Frame Sample](./images/sample_frame.png) <!-- Replace with actual image -->

---

## ⚙️ Setup & Run

```bash
# Open MATLAB
# Then in command window:
>> run code1_population_matrix.m
>> run code2_heat_equation.m
>> run code3_predator_prey_simulation.m

##📌 Disclaimer on Real-World Application
###These scripts are educational prototypes. For real ecological modeling, additional complexity such as:

-stochasticity,

-external interventions (e.g., pesticide schedules),

-habitat heterogeneity, and

-empirical calibration

-should be considered.

##📚 References
-Caswell, H. (2001). Matrix Population Models. Sinauer Associates.

-Smith, G. D. (1985). Numerical Solution of Partial Differential Equations. Oxford University Press.

-Couzin, I. D., Krause, J., et al. (2002). "Collective memory and spatial sorting in animal groups." Journal of Theoretical Biology.
