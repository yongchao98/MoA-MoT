import numpy as np

def simulate_melting(temperature, tm, cooperativity=1.0):
    """
    Simulates a simple two-state melting curve using a sigmoid function.
    
    Args:
        temperature (np.array): Array of temperature points.
        tm (float): The melting temperature (where 50% is unfolded).
        cooperativity (float): A factor that determines the sharpness of the transition.
        
    Returns:
        np.array: The fraction of unfolded molecules at each temperature.
    """
    return 1 / (1 + np.exp(cooperativity * (tm - temperature)))

# Define a temperature range for our experiment
temperatures = np.arange(40, 91, 1) # From 40C to 90C

# --- Simulation of Heterogeneity ---
# Imagine our sample has two distinct populations of molecules
# Population 1 has a melting temperature (Tm) of 60 C
tm1 = 60.0
population1_melted_fraction = simulate_melting(temperatures, tm1)

# Population 2 is more stable, with a Tm of 75 C
tm2 = 75.0
population2_melted_fraction = simulate_melting(temperatures, tm2)

# --- The "Bulk" Measurement ---
# A bulk experiment can't see the two populations. It only measures the average.
# Let's assume an equal mix (50/50).
bulk_melted_fraction = (population1_melted_fraction + population2_melted_fraction) / 2.0

# --- Print the results to show the loss of information ---
print("This simulation demonstrates how bulk measurements average out heterogeneity.")
print("We have two molecule populations with different melting points (Tm1=60C, Tm2=75C).\n")
print("A bulk experiment would observe only the average curve, obscuring the two distinct transitions.\n")

# Let's check the values at a specific temperature, for example, 62 C
temp_of_interest = 62
idx = np.where(temperatures == temp_of_interest)[0][0]

frac1 = population1_melted_fraction[idx]
frac2 = population2_melted_fraction[idx]
bulk_frac = bulk_melted_fraction[idx]

print(f"--- Analysis at Temperature = {temp_of_interest} C ---")
print(f"Individual behavior (unseen in bulk):")
print(f"  - Population 1 (Tm=60C) is {frac1:.2f} melted.")
print(f"  - Population 2 (Tm=75C) is {frac2:.2f} melted.")
print("\nBulk measurement (what is actually observed):")
print(f"The observed melted fraction is the average.")
# The following print statement shows the "equation" as requested
print(f"Equation: ({frac1:.2f} + {frac2:.2f}) / 2 = {bulk_frac:.2f}")
print(f"An experimenter would measure a melted fraction of {bulk_frac:.2f} and would not know about the two separate populations.")
