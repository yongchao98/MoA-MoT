import numpy as np

def simulate_heat_absorption(temperatures, Tm, steepness=1.0):
    """
    Simulates a heat absorption peak during melting.
    This is the derivative of the sigmoid melting curve.
    """
    exponent = (temperatures - Tm) / steepness
    # Derivative of a sigmoid function creates a bell-like curve
    return np.exp(exponent) / (steepness * (1 + np.exp(exponent))**2)

# --- Simulation Parameters ---
# Let's imagine a sample with two different types of molecules (heterogeneity)
# We will define the melting temperature (Tm) for each population.
population_1_Tm = 60.0  # Celsius
population_2_Tm = 65.0  # Celsius
fraction_pop_1 = 0.5
fraction_pop_2 = 0.5

# Define the temperature range for our simulated experiment
temperatures = np.arange(50, 75, 1)

# --- Simulate Individual and Bulk Signals ---
# Signal from Population 1 if it were measured alone
signal_1 = simulate_heat_absorption(temperatures, population_1_Tm)
# Signal from Population 2 if it were measured alone
signal_2 = simulate_heat_absorption(temperatures, population_2_Tm)

# The "Bulk" signal is the weighted average of the individual signals
# This is what an actual bulk calorimeter would measure
bulk_signal = (signal_1 * fraction_pop_1) + (signal_2 * fraction_pop_2)

# --- Print the Results ---
# The following output demonstrates the concept.
# In a real experiment, you would only see the "Bulk Signal" column.
# You would not be able to know that it's actually composed of two
# separate signals from two populations melting at 60 and 65 degrees.

print("Demonstration of Averaging in Bulk Experiments\n")
print(f"Population 1 Parameters: Tm = {population_1_Tm}°C, Fraction = {fraction_pop_1}")
print(f"Population 2 Parameters: Tm = {population_2_Tm}°C, Fraction = {fraction_pop_2}\n")

print(f"{'Temp (°C)':<12} {'Signal Pop 1':<15} {'Signal Pop 2':<15} {'Bulk Signal':<15}")
print("-" * 55)

for i, temp in enumerate(temperatures):
    print(f"{temp:<12} {signal_1[i]:<15.4f} {signal_2[i]:<15.4f} {bulk_signal[i]:<15.4f}")

print("\nConclusion: The distinct peaks from Population 1 (around 60°C) and Population 2 (around 65°C) are merged")
print("into a single, broader peak in the 'Bulk Signal'. This loss of detail is the key limitation.")
