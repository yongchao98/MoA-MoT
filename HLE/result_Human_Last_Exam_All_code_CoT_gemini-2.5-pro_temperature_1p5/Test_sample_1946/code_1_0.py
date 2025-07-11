import numpy as np
import math

# This script calculates the time it takes for a toilet paper roll to unravel completely
# when dropped from a height, using the 4th-order Runge-Kutta (RK4) method to solve the
# equations of motion.

# Step 1: Define physical constants from the problem statement.
# Gravitational acceleration (m/s^2)
g = 9.81
# Mass of the cardboard cylinder (kg)
m_cylinder = 0.020
# Total mass of the paper (kg)
m_paper_total = 0.200
# Diameter of the inner cylinder (m)
d_cylinder = 0.04
# Radius of the inner cylinder (m)
r_cylinder = d_cylinder / 2.0
# Thickness of the paper (m)
t_paper = 0.0005
# Number of times the paper is wrapped
N_wraps = 100

# Step 2: Calculate derived properties of the roll.
# Initial outer radius of the full roll (m)
R_initial = r_cylinder + N_wraps * t_paper
# Total length of the paper (m), using a continuous model for consistency.
L_total = (math.pi / t_paper) * (R_initial**2 - r_cylinder**2)
# Cross-sectional area of the paper on the roll (m^2)
A_paper_cross_section = math.pi * (R_initial**2 - r_cylinder**2)
# Area density of the paper's cross-section (kg/m^2)
sigma_area_density = m_paper_total / A_paper_cross_section
# Moment of inertia of the cardboard cylinder (approximated as a thin ring, I = mr^2)
I_cylinder = m_cylinder * r_cylinder**2

# Step 3: Define the function that represents the system of ODEs.
# This function calculates the derivatives [dy/dt, dv/dt] for a given state [y, v].
def get_derivatives(state):
    """
    Calculates the time derivatives for the ODE system.
    Input `state` is a numpy array [y, v], where y is unrolled length and v is velocity.
    Returns a numpy array [dy/dt, dv/dt].
    """
    y, v = state

    # If the roll is already unrolled, its acceleration is zero.
    if y >= L_total:
        return np.array([v, 0.0])

    # Calculate the current squared radius R(y)^2 based on unrolled length y.
    R_squared = R_initial**2 - y * t_paper / math.pi
    # Clamp the radius to the cylinder radius to avoid numerical errors near the end.
    if R_squared < r_cylinder**2:
        R_squared = r_cylinder**2
    
    # Calculate the current mass M(y).
    m_paper_remaining = sigma_area_density * math.pi * (R_squared - r_cylinder**2)
    M_current = m_cylinder + m_paper_remaining
    
    # Calculate the current moment of inertia I(y).
    I_paper_remaining = 0.5 * m_paper_remaining * (R_squared + r_cylinder**2)
    I_current = I_cylinder + I_paper_remaining
    
    # Calculate the current acceleration a(y).
    denominator = M_current + I_current / R_squared
    acceleration = (M_current * g) / denominator if denominator else 0

    return np.array([v, acceleration])

# Step 4: Run the RK4 simulation.
# Initial conditions
time = 0.0
# State vector [position y, velocity v]
state = np.array([0.0, 0.0])
# Time step for the simulation (s). A small value provides better accuracy.
dt = 0.001

# Loop until the unrolled length `y` exceeds the total length `L_total`.
while state[0] < L_total:
    # RK4 steps
    k1 = get_derivatives(state)
    k2 = get_derivatives(state + 0.5 * dt * k1)
    k3 = get_derivatives(state + 0.5 * dt * k2)
    k4 = get_derivatives(state + dt * k3)
    
    # Update the state and time
    state += (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    time += dt

# Print the final results as requested.
# The "final equation" is the condition we solved for: finding the time 't'
# when the unrolled length y(t) equals the total paper length L.
print(f"Total length of toilet paper (L): {L_total:.2f} m")
print("Final Equation: Find t such that y(t) = L")
print(f"Time = {time:.2f} s")