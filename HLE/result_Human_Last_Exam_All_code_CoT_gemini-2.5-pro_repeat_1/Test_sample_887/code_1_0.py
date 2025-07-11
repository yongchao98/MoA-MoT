# This script provides an illustrative model for the effect of bond-specific
# vibrational excitation on the reaction of F with CHD3.
# NOTE: The rate constants and factors are chosen for demonstration purposes
# to explain the chemical concept, not from precise experimental data.

def model_reaction_kinetics(k_H, k_D, num_H_bonds=1, num_D_bonds=3):
    """
    Calculates total rate and branching ratio for H abstraction.
    """
    rate_H_channel = num_H_bonds * k_H
    rate_D_channel = num_D_bonds * k_D
    total_rate = rate_H_channel + rate_D_channel
    
    # Avoid division by zero if rates are zero
    if total_rate == 0:
        branching_ratio_H = 0
    else:
        branching_ratio_H = rate_H_channel / total_rate
        
    return total_rate, branching_ratio_H

# --- Case 1: Ground Vibrational State ---
# In the ground state, C-H bonds react slightly faster than C-D bonds
# due to the kinetic isotope effect (lower zero-point energy of C-D).
# Let's assign a baseline reactivity (rate constant).
k_H_ground = 1.0  # Relative rate constant for H-abstraction per bond
k_D_ground = 0.5  # Relative rate constant for D-abstraction per bond

rate_ground, ratio_H_ground = model_reaction_kinetics(k_H_ground, k_D_ground)

print("--- Scenario 1: Ground State CHD3 Reacting with F ---")
print(f"Assumed rate constant for C-H cleavage (k_H): {k_H_ground}")
print(f"Assumed rate constant for C-D cleavage (k_D): {k_D_ground}")
print("Total Rate = (1 * k_H) + (3 * k_D)")
print(f"Total Rate = (1 * {k_H_ground}) + (3 * {k_D_ground}) = {rate_ground:.1f}")
print(f"Likelihood of H atom removal: {ratio_H_ground:.2%}")
print("-" * 55)

# --- Case 2: C-H Bond Vibrationally Excited ---
# An IR laser excites the C-H stretch. This added energy is localized
# and dramatically increases the reactivity of the C-H bond.
# The reactivity of the "cold" C-D bonds is assumed to be unaffected.
excitation_enhancement_factor = 10.0
k_H_excited = k_H_ground * excitation_enhancement_factor
k_D_excited = k_D_ground # Unchanged

rate_excited, ratio_H_excited = model_reaction_kinetics(k_H_excited, k_D_excited)
acceleration_factor = rate_excited / rate_ground

print("--- Scenario 2: C-H Excited CHD3 Reacting with F ---")
print(f"Enhanced rate constant for C-H cleavage (k_H_excited): {k_H_excited:.1f}")
print(f"Unchanged rate constant for C-D cleavage (k_D_excited): {k_D_excited}")
print("Total Rate = (1 * k_H_excited) + (3 * k_D_excited)")
print(f"Total Rate = (1 * {k_H_excited:.1f}) + (3 * {k_D_excited}) = {rate_excited:.1f}")
print(f"Likelihood of H atom removal: {ratio_H_excited:.2%}")
print("-" * 55)

print("\n--- Conclusion from Model ---")
print(f"1. The overall reaction is accelerated by a factor of {acceleration_factor:.2f}.")
print(f"2. The likelihood of H atom removal is enhanced, increasing from {ratio_H_ground:.2%} to {ratio_H_excited:.2%}.")
print("This demonstrates that exciting the C-H bond accelerates the reaction by enhancing the likelihood of H atom removal over D atoms.")
