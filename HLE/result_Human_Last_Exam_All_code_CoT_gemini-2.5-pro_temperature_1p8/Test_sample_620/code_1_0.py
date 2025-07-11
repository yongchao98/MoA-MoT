import math

def calculate_active_dimer(total_enzyme, kd):
    """
    Calculates the concentration of active dimer based on total enzyme concentration and Kd.
    
    Equilibrium: Dimer <-> 2 Monomer
    Dissociation constant: Kd = [Monomer]^2 / [Dimer]
    Total enzyme (monomer units): E_total = [Monomer] + 2 * [Dimer]
    
    This solves the quadratic equation: 2*D + sqrt(Kd*D) - E_total = 0
    
    Args:
        total_enzyme (float): The total concentration of enzyme monomer units.
        kd (float): The dissociation constant for the dimer.
        
    Returns:
        tuple: A tuple containing (concentration of dimer, fraction of enzyme in dimer form).
    """
    # Using the quadratic formula to solve for sqrt([Dimer])
    # x = (-b + sqrt(b^2 - 4ac)) / 2a
    # where x = sqrt([Dimer]), a = 2, b = sqrt(Kd), c = -E_total
    sqrt_dimer = (-math.sqrt(kd) + math.sqrt(kd + 8 * total_enzyme)) / 4
    dimer_concentration = sqrt_dimer**2
    
    # Active fraction is the concentration of monomers in the dimer form
    # divided by the total concentration of monomers.
    active_fraction = (2 * dimer_concentration) / total_enzyme
    
    return dimer_concentration, active_fraction

# --- Simulation Parameters ---
# We assume a dissociation constant (Kd). A higher Kd means more dissociation.
# Chilling the enzyme likely increases its Kd. Let's assume Kd = 10 (in arbitrary units).
KD_VALUE = 10.0

# Case 1: Low enzyme concentration
low_enzyme_concentration = 1.0

# Case 2: High enzyme concentration (100-fold higher)
high_enzyme_concentration = 100.0

# --- Calculations ---
dimer_low, fraction_low = calculate_active_dimer(low_enzyme_concentration, KD_VALUE)
dimer_high, fraction_high = calculate_active_dimer(high_enzyme_concentration, KD_VALUE)

# --- Output the results ---
print("--- Modeling Enzyme Dimerization ---")
print(f"Assumed Dissociation Constant (Kd): {KD_VALUE}\n")
print("This model shows how increasing total enzyme concentration pushes the equilibrium towards the active dimer form.\n")

print(f"Scenario 1: Low Enzyme Concentration")
print(f"Total Enzyme Concentration = {low_enzyme_concentration}")
print(f"Calculated Active Dimer Concentration = {dimer_low:.4f}")
print(f"Fraction of Enzyme in Active Dimer Form = {fraction_low:.2%}\n")

print(f"Scenario 2: High Enzyme Concentration")
print(f"Total Enzyme Concentration = {high_enzyme_concentration}")
print(f"Calculated Active Dimer Concentration = {dimer_high:.4f}")
print(f"Fraction of Enzyme in Active Dimer Form = {fraction_high:.2%}\n")

print("Conclusion: Increasing enzyme concentration from 1.0 to 100.0 increased the")
print("fraction of active enzyme from ~15% to ~80%. This would solve the lag phase issue.")
print("Therefore, the best troubleshooting step is C.")
