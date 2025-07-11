# Step 1: Define the known values from the problem statement.
mass_initial_solution = 10.0  # g
mass_fraction_initial_salt = 0.10
plate_mass_decrease = 0.172  # g
# A is divalent
mass_fraction_final_salt = 0.1152

# Step 2: Calculate initial and final masses based on the given data.
mass_initial_salt = mass_initial_solution * mass_fraction_initial_salt
mass_final_solution = mass_initial_solution + plate_mass_decrease
# As shown in the thinking process, we use the value for the final salt mass
# that is stoichiometrically consistent with the other mass data.
mass_final_salt = mass_initial_salt + plate_mass_decrease

# Step 3: Establish the relationship between the atomic weights of Metal A (Ar_A) and Metal M (Ar_M).
# We derived the formula: Ar_A = 1.172 * Ar_M + 12.212
# Let's test a common divalent metal for M, like Copper (Cu).
Ar_M_copper = 63.546  # Atomic weight of Copper

# Calculate the theoretical Ar_A based on M being Copper
Ar_A_calculated = 1.172 * Ar_M_copper + 12.212
Ar_Sr = 87.62 # Atomic weight of Strontium

# Step 4: Final determination and outputting the equation.
# The calculated Ar_A is very close to that of Strontium (Sr), a divalent metal
# which is more reactive than Copper. This confirms our identification.
metal_A_name = "Strontium (Sr)"
metal_M_name = "Copper (Cu)"
metal_A_symbol = "Sr"
unknown_chloride_formula = "CuCl2"
final_salt_formula = "SrCl2"
metal_M_symbol = "Cu"

# Print the determined metals and the final reaction equation.
print(f"The analysis suggests that Metal A is {metal_A_name} and the unknown chloride was {metal_M_name} chloride.")
print(f"The calculated atomic weight for Metal A was {Ar_A_calculated:.2f}, which is very close to the actual atomic weight of Strontium ({Ar_Sr}).")
print("\nThe balanced chemical equation for the reaction is:")

# Output the final equation, including each number/symbol as requested
print(f"{metal_A_symbol} + {unknown_chloride_formula} -> {final_salt_formula} + {metal_M_symbol}")

# Output final answer in the specified format
# The question asks to "Determine the metal and write the equation".
# The primary answer is the metal A.
<<<Strontium>>>