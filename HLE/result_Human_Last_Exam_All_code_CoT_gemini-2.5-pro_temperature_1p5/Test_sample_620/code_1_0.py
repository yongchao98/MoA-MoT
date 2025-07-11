import numpy as np

def simulate_enzyme_kinetics(enzyme_concentration, substrate_concentration, Km, kcat, duration_s, time_step_s):
    """
    Simulates product formation in an enzyme-catalyzed reaction over time.

    The function models the change in substrate and product concentration based on the
    Michaelis-Menten equation: v = (kcat * [E] * [S]) / (Km + [S])

    Args:
        enzyme_concentration (float): The initial concentration of the enzyme [E].
        substrate_concentration (float): The initial concentration of the substrate [S].
        Km (float): The Michaelis constant.
        kcat (float): The turnover number (catalytic rate constant).
        duration_s (int): The total time of the simulation in seconds.
        time_step_s (int): The interval between measurements in seconds.

    Returns:
        A list of tuples, where each tuple contains (time, product_concentration).
    """
    time_points = np.arange(0, duration_s + time_step_s, time_step_s)
    results = []

    current_S = substrate_concentration
    current_P = 0.0

    # Ensure we don't try to get a value for a time point that isn't calculated
    # This loop calculates the state at each discrete time step.
    calculated_P = {0.0: 0.0}
    for t in np.arange(0, duration_s, 1): # calculate every second
        if current_S > 0:
            rate = (kcat * enzyme_concentration * current_S) / (Km + current_S)
            delta_S = rate * 1 # change over 1 second
            current_S -= delta_S
            current_P += delta_S
        calculated_P[t+1] = current_P


    for t in time_points:
        results.append((t, calculated_P[t]))

    return results

# --- Shared Parameters for the Assay ---
KM_VALUE = 50.0   # Michaelis constant in uM
KCAT_VALUE = 20.0 # Turnover number in 1/s
S0_VALUE = 100.0  # Initial substrate concentration in uM
ASSAY_TIME = 60 # seconds

# --- Scenario 1: The Problem (High Enzyme Concentration) ---
# The non-linear curve suggests the enzyme concentration is too high,
# causing rapid substrate depletion.
high_enzyme_conc = 1.0 # uM

print("--- Problem Scenario: High Enzyme Concentration ---")
print(f"Simulating with [E] = {high_enzyme_conc} uM, [S] = {S0_VALUE} uM, Km = {KM_VALUE} uM, kcat = {KCAT_VALUE} 1/s")

high_enzyme_results = simulate_enzyme_kinetics(high_enzyme_conc, S0_VALUE, KM_VALUE, KCAT_VALUE, ASSAY_TIME, 5)

print("\nProduct vs. Time (High [E])")
print("Time (s) | Product (uM) | Rate (uM/s, avg over interval)")
print("---------------------------------------------------------")
last_product = 0.0
last_time = 0
for time, product in high_enzyme_results:
    if time > last_time:
        rate = (product - last_product) / (time - last_time)
        print(f"{time:<8} | {product:<12.2f} | {rate:<.2f}")
    else:
        print(f"{time:<8} | {product:<12.2f} | N/A")
    last_product = product
    last_time = time
print("\nObservation: The rate drops significantly, showing the curve is not linear.\n")


# --- Scenario 2: The Solution (Decrease Enzyme Concentration) ---
# The recommended troubleshooting step is to decrease the enzyme concentration.
low_enzyme_conc = 0.1 # uM (10-fold dilution)

print("--- Solution: Decreased Enzyme Concentration ---")
print(f"Simulating with [E] = {low_enzyme_conc} uM, [S] = {S0_VALUE} uM, Km = {KM_VALUE} uM, kcat = {KCAT_VALUE} 1/s")

low_enzyme_results = simulate_enzyme_kinetics(low_enzyme_conc, S0_VALUE, KM_VALUE, KCAT_VALUE, ASSAY_TIME, 5)

print("\nProduct vs. Time (Low [E])")
print("Time (s) | Product (uM) | Rate (uM/s, avg over interval)")
print("---------------------------------------------------------")
last_product = 0.0
last_time = 0
for time, product in low_enzyme_results:
    if time > last_time:
        rate = (product - last_product) / (time - last_time)
        print(f"{time:<8} | {product:<12.2f} | {rate:<.2f}")
    else:
        print(f"{time:<8} | {product:<12.2f} | N/A")
    last_product = product
    last_time = time
print("\nObservation: The rate is much more constant, indicating a linear phase suitable for analysis.")