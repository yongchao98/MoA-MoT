import numpy as np

def simulate_kinetics(E0, S0, Km, k_cat, time_points):
    """
    Simulates product formation over time using numerical integration
    of the Michaelis-Menten equation.
    V = (Vmax * [S]) / (Km + [S]), where Vmax = k_cat * [E]
    """
    Vmax = k_cat * E0
    S = S0
    P = 0
    # Use a smaller time step for more accurate numerical simulation
    dt = 0.1
    total_time = time_points[-1]
    
    product_over_time = []
    current_time_point_index = 0

    # Simulate in small steps (dt) but only record at specified time_points
    for step in range(int(total_time / dt) + 1):
        t = step * dt
        if current_time_point_index < len(time_points) and t >= time_points[current_time_point_index]:
            product_over_time.append(P)
            current_time_point_index += 1

        rate = (Vmax * S) / (Km + S)
        delta_P = rate * dt
        P += delta_P
        S -= delta_P
        if S < 0:
            S = 0
            
    return np.array(product_over_time)

# --- Main Analysis ---

# Plan: The problem describes a scenario where the Product vs. Time plot is a curve
# with no linear phase. This typically occurs when the enzyme concentration is too high,
# leading to rapid consumption of the substrate. The best troubleshooting step is to
# slow the reaction down to make the initial, linear phase measurable. The most direct
# way to do this is by reducing the enzyme concentration.
#
# We will simulate two scenarios to demonstrate this:
# 1. A "Problematic Assay" with a high enzyme concentration.
# 2. A "Corrected Assay" with a 50-fold lower enzyme concentration.

print("--- Troubleshooting an Enzyme Kinetics Assay ---")
print("Problem: Product vs. Time plot is curved, not linear, indicating the reaction is too fast.\n")

# --- Assay Parameters ---
S0 = 100.0  # Initial Substrate Concentration (µM)
Km = 20.0   # Michaelis Constant (µM)
k_cat = 50.0 # Turnover Number (1/s)
# Time points for measurement (e.g., every 5 seconds for 1 minute)
time_points = np.arange(0, 61, 5) 

# --- Scenario 1: Problematic Assay (High Enzyme Concentration) ---
E0_high = 1.0 # µM
Vmax_high = k_cat * E0_high
product_high_E = simulate_kinetics(E0_high, S0, Km, k_cat, time_points)

print("--- SCENARIO 1: Problematic Assay (High [Enzyme]) ---")
print(f"Initial Substrate [S] = {S0} µM")
print(f"Enzyme Concentration [E] = {E0_high} µM")
print(f"Resulting Vmax = {Vmax_high} µM/s")
print("\nProduct Formation (µM):")
print("Time (s) | Product (µM)")
print("-----------------------")
for t, p in zip(time_points, product_high_E):
    print(f"  {t:2d}     |   {p:5.1f}")
print("\nObservation: The reaction is nearly complete by 15 seconds. The rate drops very fast, so there is no linear phase to measure.\n")

# --- Scenario 2: Corrected Assay (Decreased Enzyme Concentration) ---
# Decrease the enzyme concentration by 50-fold to slow the reaction
E0_low = E0_high / 50.0 
Vmax_low = k_cat * E0_low
product_low_E = simulate_kinetics(E0_low, S0, Km, k_cat, time_points)

print("--- SCENARIO 2: Corrected Assay (Low [Enzyme]) ---")
print(f"Action: Decrease Enzyme Concentration from {E0_high} µM to {E0_low} µM.")
print(f"New Enzyme Concentration [E] = {E0_low} µM")
print(f"Resulting Vmax = {Vmax_low} µM/s")
print("\nProduct Formation (µM):")
print("Time (s) | Product (µM) | Rate in Interval (µM/s)")
print("-------------------------------------------------")
# Calculate rate between points to show linearity
rates_low_E = np.diff(product_low_E) / np.diff(time_points)
print(f"  {time_points[0]:2d}     |   {product_low_E[0]:5.1f}    |")
for i in range(len(rates_low_E)):
    print(f"  {time_points[i+1]:2d}     |   {product_low_E[i+1]:5.1f}    |         {rates_low_E[i]:.2f}")

print("\nObservation: Product formation is now slower and clearly linear over 60 seconds.")
print("The rate is constant in the early time points, allowing for accurate initial velocity measurement.")
print("\nConclusion: Decreasing the enzyme concentration is the correct troubleshooting step.")