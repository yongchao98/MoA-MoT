import math

def solve_dimer_equilibrium(total_enzyme_conc, Ka):
    """
    Calculates the concentration of monomer and active dimer at equilibrium.

    The equilibrium is 2*M <-> D, where M is monomer and D is dimer.
    The total enzyme concentration in terms of monomers is [Et] = [M] + 2*[D].
    The association constant is Ka = [D] / [M]^2.

    Substituting [D] gives a quadratic equation for [M]:
    2*Ka*[M]^2 + [M] - [Et] = 0
    
    This function solves this equation for [M] and then calculates [D].
    """
    
    # Coefficients for the quadratic formula: a*x^2 + b*x + c = 0
    a = 2 * Ka
    b = 1
    c = -total_enzyme_conc

    # Use the quadratic formula to solve for [Monomer]
    # x = (-b + sqrt(b^2 - 4ac)) / 2a
    # We use the positive root because concentration cannot be negative.
    discriminant = b**2 - 4*a*c
    if discriminant < 0:
        return "Error: No real solution exists. Check parameters."
        
    monomer_conc = (-b + math.sqrt(discriminant)) / (2*a)
    
    # Calculate dimer concentration from the monomer concentration
    dimer_conc = Ka * monomer_conc**2
    
    # Calculate the percentage of the enzyme that is in the active dimer form
    active_dimer_fraction = (2 * dimer_conc) / total_enzyme_conc * 100

    print(f"For a total enzyme concentration of {total_enzyme_conc:.2e} M:")
    print(f"  The equation to solve is: {a:.2e}*[M]^2 + {b}*[M] - {total_enzyme_conc:.2e} = 0")
    print(f"  Calculated [Monomer] = {monomer_conc:.2e} M")
    print(f"  Calculated [Dimer] = {dimer_conc:.2e} M")
    print(f"  This means {active_dimer_fraction:.1f}% of the enzyme is in the active dimer form.")
    print("-" * 50)

# --- Simulation ---

# Assume an association constant (Ka) for the dimer. 
# A lower Ka means the dimer is less stable and more likely to dissociate.
# Let's assume the chilling on ice has lowered the effective Ka.
# Ka units are M^-1
association_constant = 1e6 # 1 / uM

# Case 1: Original (Low) Enzyme Concentration
low_enzyme_conc = 2e-6 # 2 uM
print(">>> CASE 1: LOW ENZYME CONCENTRATION <<<")
solve_dimer_equilibrium(low_enzyme_conc, association_constant)

# Case 2: Increased Enzyme Concentration (Choice C)
# Let's increase the concentration by 5-fold
high_enzyme_conc = 10e-6 # 10 uM
print(">>> CASE 2: INCREASED ENZYME CONCENTRATION <<<")
solve_dimer_equilibrium(high_enzyme_conc, association_constant)

print("\nCONCLUSION:")
print("Increasing the total enzyme concentration significantly increases the percentage")
print("of the enzyme in the active, dimeric form. This would reduce or eliminate")
print("the lag phase in the kinetics assay, leading to a more linear plot.")
print("Therefore, increasing the enzyme concentration is the correct action.")

<<<C>>>