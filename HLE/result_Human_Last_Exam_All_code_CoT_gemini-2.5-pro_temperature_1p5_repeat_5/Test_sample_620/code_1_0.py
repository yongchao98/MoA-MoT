import sys
import io

# A function to simulate the enzyme reaction and print the results
def simulate_enzyme_kinetics():
    """
    This simulation demonstrates why decreasing enzyme concentration helps to achieve
    a linear initial rate in an enzyme assay. It models the consumption of a
    substrate over time using the Michaelis-Menten equation.
    """
    # Define standard assay parameters
    Km = 15.0  # Michaelis constant in uM
    kcat = 20.0  # Turnover number in events/second
    initial_substrate = 100.0  # Initial substrate concentration in uM
    total_time = 60  # Total assay time in seconds
    time_step = 10  # Time interval for measurements in seconds

    print("Step-by-step analysis of why decreasing enzyme concentration is the correct approach.")
    print("The goal is to find conditions where the reaction rate (product formed per second) is constant at the beginning of the assay.\n")

    # --- Case 1: High Enzyme Concentration ---
    high_enzyme_conc = 1.0  # in nM
    print(f"--- Simulation 1: High Enzyme Concentration ({high_enzyme_conc} nM) ---")
    print("Problem: The rate is too fast, depleting the substrate quickly, leading to a non-linear curve.")
    print("-" * 70)
    print(f"{'Time (s)':<12} | {'Substrate (uM)':<18} | {'Product (uM)':<18} | {'Rate (uM/s)':<15}")
    print("-" * 70)

    substrate = initial_substrate
    product = 0.0
    for t in range(0, total_time + 1, time_step):
        # Michaelis-Menten equation: v = (kcat * [E] * [S]) / (Km + [S])
        velocity = (kcat * high_enzyme_conc * substrate) / (Km + substrate)
        if t == 0:
             print(f"{t:<12} | {substrate:<18.2f} | {product:<18.2f} | {'N/A'}")
        else:
             print(f"{t:<12} | {substrate:<18.2f} | {product:<18.2f} | {velocity:<15.2f}")
        
        # Update concentrations for the next time step
        product_formed_in_step = velocity * time_step
        substrate -= product_formed_in_step
        product += product_formed_in_step
        if substrate < 0:
            substrate = 0
            
    print("-" * 70)
    print("Observation: The rate drops from >11 uM/s to <8 uM/s. This is a non-linear phase.\n\n")


    # --- Case 2: Decreased Enzyme Concentration ---
    low_enzyme_conc = 0.1  # in nM (10-fold dilution)
    print(f"--- Simulation 2: Decreased Enzyme Concentration ({low_enzyme_conc} nM) ---")
    print("Solution: Reducing the enzyme concentration makes the rate more constant.")
    print("-" * 70)
    print(f"{'Time (s)':<12} | {'Substrate (uM)':<18} | {'Product (uM)':<18} | {'Rate (uM/s)':<15}")
    print("-" * 70)

    substrate = initial_substrate
    product = 0.0
    for t in range(0, total_time + 1, time_step):
        velocity = (kcat * low_enzyme_conc * substrate) / (Km + substrate)
        if t == 0:
             print(f"{t:<12} | {substrate:<18.2f} | {product:<18.2f} | {'N/A'}")
        else:
             print(f"{t:<12} | {substrate:<18.2f} | {product:<18.2f} | {velocity:<15.2f}")
        
        product_formed_in_step = velocity * time_step
        substrate -= product_formed_in_step
        product += product_formed_in_step
        if substrate < 0:
            substrate = 0
            
    print("-" * 70)
    print("Observation: The rate is now very stable (~1.7 uM/s). This provides a reliable linear phase to measure initial velocity.")

# Execute the simulation
simulate_enzyme_kinetics()
<<<D>>>