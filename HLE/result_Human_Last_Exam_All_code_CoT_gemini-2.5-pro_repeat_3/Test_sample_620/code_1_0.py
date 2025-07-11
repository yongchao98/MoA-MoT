import pandas as pd

def simulate_reaction(enzyme_conc, initial_substrate, Km, kcat, time_points):
    """
    Simulates product formation over time using the Michaelis-Menten equation,
    accounting for substrate depletion.
    """
    product_formed_list = []
    substrate_conc = initial_substrate
    
    # Calculate product at each time point
    for t in time_points:
        if t == 0:
            product_formed = 0.0
        else:
            # For this simulation, we calculate the change over each 1-second interval (dt=1)
            # Velocity v = (kcat * [E] * [S]) / (Km + [S])
            velocity = (kcat * enzyme_conc * substrate_conc) / (Km + substrate_conc)
            substrate_consumed = velocity * 1.0  # dt = 1s
            substrate_conc -= substrate_consumed
            if substrate_conc < 0:
                substrate_conc = 0
            
            product_formed = initial_substrate - substrate_conc
            
        product_formed_list.append(round(product_formed, 2))
        
    return product_formed_list

def main():
    # --- Enzyme and Assay Parameters ---
    Km = 10.0  # Michaelis constant (in uM)
    kcat = 50.0  # Turnover number (in events/sec)
    initial_substrate = 100.0  # Initial substrate concentration (in uM)
    time_seconds = list(range(11)) # 0 to 10 seconds
    
    # Scenario 1: High Enzyme Concentration
    high_enzyme_conc = 1.0  # in uM
    
    # Scenario 2: Decreased Enzyme Concentration
    low_enzyme_conc = 0.1  # in uM (10-fold decrease)

    print("--- Problem Analysis ---")
    print("A non-linear Product vs. Time plot often means the reaction is too fast, depleting the substrate.")
    print("This simulation will show why decreasing enzyme concentration helps linearize the plot.\n")

    print("--- Simulation Parameters ---")
    print(f"Km = {Km} uM")
    print(f"kcat = {kcat}/sec")
    print(f"Initial [Substrate] = {initial_substrate} uM\n")

    # --- Run Simulations ---
    product_high_E = simulate_reaction(high_enzyme_conc, initial_substrate, Km, kcat, time_seconds)
    product_low_E = simulate_reaction(low_enzyme_conc, initial_substrate, Km, kcat, time_seconds)
    
    # --- Display Results ---
    results_df = pd.DataFrame({
        'Time (s)': time_seconds,
        'Product (High [E])': product_high_E,
        'Product (Low [E])': product_low_E
    })
    
    print("--- Simulated Product Formation (uM) ---")
    print(results_df.to_string(index=False))
    print("\n--- Analysis of Initial Rates ---")
    
    # High [E] analysis
    print("\nScenario 1: High Enzyme Concentration")
    print(f"Enzyme Concentration [E] = {high_enzyme_conc} uM")
    rate_high_E_t1 = product_high_E[1] / 1.0
    # Show calculation for the first second
    print(f"Initial Velocity Equation: v = (kcat * [E] * [S]) / (Km + [S])")
    print(f"v (at t=0) = ({kcat} * {high_enzyme_conc} * {initial_substrate}) / ({Km} + {initial_substrate}) = {rate_high_E_t1:.2f} uM/s")
    
    rate_high_E_t5 = (product_high_E[5] - product_high_E[4]) / 1.0
    print(f"Rate between t=4s and t=5s: {rate_high_E_t5:.2f} uM/s")
    linearity_check_high = (rate_high_E_t5 / rate_high_E_t1) * 100
    print(f"Result: The rate drops significantly, showing the plot is a curve. The rate at 5s is only {linearity_check_high:.1f}% of the initial rate.")
    
    # Low [E] analysis
    print("\nScenario 2: Decreased Enzyme Concentration")
    print(f"Enzyme Concentration [E] = {low_enzyme_conc} uM")
    rate_low_E_t1 = product_low_E[1] / 1.0
    # Show calculation for the first second
    print(f"Initial Velocity Equation: v = (kcat * [E] * [S]) / (Km + [S])")
    print(f"v (at t=0) = ({kcat} * {low_enzyme_conc} * {initial_substrate}) / ({Km} + {initial_substrate}) = {rate_low_E_t1:.2f} uM/s")

    rate_low_E_t5 = (product_low_E[5] - product_low_E[4]) / 1.0
    print(f"Rate between t=4s and t=5s: {rate_low_E_t5:.2f} uM/s")
    linearity_check_low = (rate_low_E_t5 / rate_low_E_t1) * 100
    print(f"Result: The rate is much more stable. The rate at 5s is {linearity_check_low:.1f}% of the initial rate, indicating a more linear plot.")

if __name__ == '__main__':
    main()