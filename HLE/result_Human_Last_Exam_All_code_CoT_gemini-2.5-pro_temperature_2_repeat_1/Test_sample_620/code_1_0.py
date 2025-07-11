import numpy as np
import matplotlib.pyplot as plt

def simulate_enzyme_kinetics(E_total, S_initial, Km, kcat, time_points):
    """
    Simulates product formation over time using the Michaelis-Menten model.
    Vmax is a function of total enzyme concentration (E_total) and kcat.
    """
    S = S_initial
    P = 0
    product_over_time = []
    
    dt = time_points[1] - time_points[0] # time step
    Vmax = kcat * E_total
    
    for t in time_points:
        # Michaelis-Menten equation for reaction rate (velocity)
        v = (Vmax * S) / (Km + S)
        
        # Amount of substrate consumed in this time step
        delta_S = v * dt
        
        # Update substrate and product concentrations
        # Ensure substrate doesn't go below zero
        if S >= delta_S:
            S -= delta_S
            P += delta_S
        else:
            P += S
            S = 0
            
        product_over_time.append(P)
        
    return time_points, product_over_time

# --- Simulation Parameters ---
S_initial = 100  # Initial Substrate Concentration (uM)
Km = 20          # Michaelis constant (uM)
kcat = 10        # Turnover number (per second)
time = np.linspace(0, 60, 300) # 60 seconds of reaction time

# --- Scenario 1: High Enzyme Concentration ---
E_high = 2.0  # uM
print("Scenario 1: High Enzyme Concentration")
print("This scenario simulates the user's problem.")
print("The enzyme consumes the substrate so quickly that the plot of Product vs. Time is a curve, not a line.")
print("It is difficult to measure a reliable initial rate.\n")
t1, p1 = simulate_enzyme_kinetics(E_high, S_initial, Km, kcat, time)

# --- Scenario 2: Decreased Enzyme Concentration ---
E_low = 0.2   # uM (10-fold less)
print("Scenario 2: Decreased Enzyme Concentration (The Solution)")
print("By decreasing the enzyme concentration, the overall reaction slows down.")
print("This results in a clear linear phase at the beginning of the reaction, suitable for measuring the initial velocity.")
t2, p2 = simulate_enzyme_kinetics(E_low, S_initial, Km, kcat, time)

# --- Plotting the results ---
plt.figure(figsize=(10, 6))
plt.plot(t1, p1, label=f'High [E] = {E_high} uM', color='red')
plt.plot(t2, p2, label=f'Decreased [E] = {E_low} uM', color='blue')

# Highlight the initial linear phase in the good assay
plt.plot(t2[:50], p2[:50], color='green', linewidth=4, label='Ideal Linear Phase')

plt.title('Troubleshooting Non-Linear Enzyme Kinetics')
plt.xlabel('Time (seconds)')
plt.ylabel('Product Formed (uM)')
plt.legend()
plt.grid(True)
print("\nPlotting the simulation. The plot window needs to be closed to end the script.")
plt.show()
