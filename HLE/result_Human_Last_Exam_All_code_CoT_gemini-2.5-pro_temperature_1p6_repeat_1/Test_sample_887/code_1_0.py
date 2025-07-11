import numpy as np

# This script models the effect of vibrational excitation on chemical reactivity.

# --- Step 1: Define Base Reactivities ---
# In an unexcited molecule, the C-H bond is slightly more reactive than a C-D bond
# due to the kinetic isotope effect (lower zero-point energy).
# Let's assign some arbitrary reactivity values.
base_reactivity_H = 1.5
base_reactivity_D = 1.0

# --- Step 2: Define Vibrational Excitation Effect ---
# Exciting the C-H bond with a laser dramatically increases its reactivity.
# Experiments show this can enhance the reaction rate by a factor of 100 or more.
vibrational_enhancement_factor = 100.0
reactivity_H_excited = base_reactivity_H * vibrational_enhancement_factor

# --- Step 3: Analyze the Unexcited Molecule (Ground State) ---
# The molecule is CHD3, so it has one C-H bond and three C-D bonds.
total_reactivity_unexcited = base_reactivity_H + (3 * base_reactivity_D)
prob_H_removal_unexcited = base_reactivity_H / total_reactivity_unexcited
prob_D_removal_unexcited = (3 * base_reactivity_D) / total_reactivity_unexcited

print("--- Analysis for Unexcited CHD3 Molecule ---")
print(f"Reactivity of C-H bond: {base_reactivity_H}")
print(f"Reactivity of each C-D bond: {base_reactivity_D}")
print(f"Total molecular reactivity (proportional to rate): {total_reactivity_unexcited:.2f}")
print(f"Likelihood of H atom removal: {prob_H_removal_unexcited:.2%}")
print(f"Likelihood of D atom removal: {prob_D_removal_unexcited:.2%}\n")


# --- Step 4: Analyze the Excited Molecule ---
# Now we excite the C-H bond. The C-D bonds remain in their ground state.
total_reactivity_excited = reactivity_H_excited + (3 * base_reactivity_D)
prob_H_removal_excited = reactivity_H_excited / total_reactivity_excited
prob_D_removal_excited = (3 * base_reactivity_D) / total_reactivity_excited

print("--- Analysis for CHD3 with Excited C-H Bond ---")
print(f"Reactivity of EXCITED C-H bond: {reactivity_H_excited}")
print(f"Reactivity of each C-D bond: {base_reactivity_D}")
print(f"Total molecular reactivity (proportional to rate): {total_reactivity_excited:.2f}")
print(f"Likelihood of H atom removal: {prob_H_removal_excited:.2%}")
print(f"Likelihood of D atom removal: {prob_D_removal_excited:.2%}\n")


# --- Step 5: Conclusion ---
rate_acceleration = total_reactivity_excited / total_reactivity_unexcited
print("--- Conclusion ---")
print(f"The overall reaction rate is accelerated by a factor of {rate_acceleration:.1f}.")
print(f"The likelihood of removing an H atom is enhanced from {prob_H_removal_unexcited:.2%} to {prob_H_removal_excited:.2%}.")
print("\nThis demonstrates that exciting the C-H bond accelerates the reaction and enhances the likelihood of H atom removal over D atoms.")