def solve_kinetics_problem():
    """
    This function simulates an enzyme kinetics assay to demonstrate why
    decreasing enzyme concentration helps to linearize the product vs. time plot.
    """
    # --- Define Enzyme and Assay Parameters ---
    # Michaelis constant (uM), representing substrate binding affinity
    KM = 50.0
    # Turnover number (per second), representing catalytic speed per enzyme molecule
    KCAT = 20.0
    # Initial substrate concentration (uM)
    S0 = 100.0

    # --- Simulation Time Parameters ---
    TOTAL_TIME = 60 # Total simulation time in seconds
    TIME_STEP = 10  # Time interval for measurements in seconds

    def simulate_reaction(enzyme_conc):
        """
        Simulates the reaction for a given enzyme concentration and prints the results.
        """
        # Vmax is the maximum reaction rate, directly proportional to enzyme concentration
        vmax = KCAT * enzyme_conc
        
        # Initialize starting conditions
        s_current = S0
        p_current = 0.0
        
        print(f"\n--- Simulating with Enzyme Concentration = {enzyme_conc:.2f} uM (Vmax = {vmax:.2f} uM/s) ---")
        print("This scenario corresponds to the rate if we use this 'Equation': Rate = (Vmax * [S]) / (Km + [S])")
        print("-" * 60)
        print(f"{'Time (s)':<12} | {'Product (uM)':<15} | {'Rate in Interval (uM/s)':<25}")
        print("-" * 60)
        
        last_p = 0.0
        # Loop through time, taking measurements at each TIME_STEP
        for t in range(0, TOTAL_TIME + 1, TIME_STEP):
            if t == 0:
                print(f"{t:<12} | {p_current:<15.2f} | {'N/A':<25}")
                continue

            # Michaelis-Menten equation to calculate the rate at the start of the interval
            rate = (vmax * s_current) / (KM + s_current)
            
            # Calculate product change over the time step
            p_delta = rate * TIME_STEP
            
            # Update concentrations for the next step
            p_current += p_delta
            s_current -= p_delta
            s_current = max(0, s_current) # Substrate cannot be negative
            
            # Calculate the average rate over the last interval to check for linearity
            interval_rate = (p_current - last_p) / TIME_STEP
            
            # Print each number for the final equation 'components' at each step
            print(f"{t:<12} | {p_current:<15.2f} | {interval_rate:<25.3f}")
            last_p = p_current

    print("### Scenario 1: High Enzyme Concentration (The Problem) ###")
    print("Notice how the 'Rate in Interval' decreases quickly, causing a non-linear plot.")
    # Represents the problematic assay with a high enzyme concentration
    e_high = 1.0 
    simulate_reaction(e_high)

    print("\n### Scenario 2: Decreased Enzyme Concentration (The Solution) ###")
    print("Notice how the 'Rate in Interval' is much more stable, indicating a linear phase.")
    # Represents the corrected assay with a decreased enzyme concentration
    e_low = 0.1
    simulate_reaction(e_low)

solve_kinetics_problem()