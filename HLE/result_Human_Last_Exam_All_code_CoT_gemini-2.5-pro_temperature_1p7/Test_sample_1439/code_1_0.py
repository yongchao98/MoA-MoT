import sympy

# --- Symbolic Calculation of the Critical Exponent ν ---

# This script will determine the order in the coupling constant 'u' at which the
# critical exponent 'ν' receives its first non-vanishing correction in φ⁴ scalar field theory.
# We use the perturbative ε-expansion framework and symbolic mathematics to derive the result.

print("Derivation of the first-order correction to the critical exponent ν:")
print("=" * 70)

# Define the symbols we will use in our calculation
u, epsilon, N = sympy.symbols('u epsilon N')

# Step 1: Define the one-loop beta function β(u) for the O(N) model
print("Step 1: The one-loop beta function β(u) for the O(N) model is:")
# Using a common normalization scheme for the coupling constant u
beta_u = u * (-epsilon + (N + 8) / 6 * u)
print(f"   β(u) = u * (-ε + (N+8)/6 * u)")
print("-" * 70)

# Step 2: Find the non-trivial fixed point u*
print("Step 2: Find the non-trivial Wilson-Fisher fixed point u* where β(u*) = 0.")
# Solve beta_u = 0 for u. The solutions are 0 (trivial) and the one we want.
fixed_point_sols = sympy.solve(beta_u, u)
u_star = fixed_point_sols[1] # The non-zero solution u*
print(f"   The non-trivial solution is u* = {u_star}")
print("   Note that the fixed point coupling u* is of the first order in ε.")
print("-" * 70)

# Step 3: State the relation for ν and the one-loop result for γ_φ²(u)
print("Step 3: State the relationship between ν, the exponent, and γ_φ², the anomalous dimension.")
print("   The exponent ν is found from the relation: 1/ν = 2 - γ_φ²(u*).")
print("   To one-loop, the anomalous dimension is linear in the coupling u:")
gamma_phi_sq_expr = (N + 2) / 3 * u
print(f"   γ_φ²(u) = {(N+2)}/3 * u")
print("-" * 70)

# Step 4: Calculate ν at the fixed point, specializing to φ⁴ theory (N=1)
print("Step 4: Calculate ν at the fixed point for φ⁴ theory (N=1).")

# First, substitute N=1 into the expressions for u* and γ_φ²(u)
u_star_n1 = u_star.subs(N, 1)
gamma_expr_n1 = gamma_phi_sq_expr.subs(N, 1)
print(f"   For N=1, u* becomes: {sympy.simplify(u_star_n1)}")

# Next, evaluate γ_φ² at the fixed point u* for N=1
gamma_at_fp = gamma_expr_n1.subs(u, u_star_n1)
print(f"   Evaluating γ_φ²(u*) at N=1 gives: {sympy.simplify(gamma_at_fp)}")

# Now, construct the expression for 1/ν
nu_inverse_expr = 2 - gamma_at_fp
print(f"   The expression for 1/ν is: 2 - {sympy.simplify(gamma_at_fp)} = {sympy.simplify(nu_inverse_expr)}")

# Finally, we find ν by inverting this and expanding for small ε
nu_expr = 1 / nu_inverse_expr
nu_series = nu_expr.series(epsilon, 0, 2).removeO() # Expand to O(ε^2) and remove the O term

# Extract the coefficients to display the final equation clearly
c0 = nu_series.coeff(epsilon, 0)
c1 = nu_series.coeff(epsilon, 1)

print("\n   The final equation for ν expanded to first order in ε is:")
print(f"   ν = {c0} + {c1}*ε")
print("-" * 70)

# Step 5: Conclusion
print("Conclusion:")
print("The zeroth-order term ν = 1/2 is the mean-field theory result.")
print("The first non-vanishing correction to this value is the term linear in ε.")
print("Since we found that the fixed point coupling u* is also linear in ε (u* ∝ ε), a correction of order O(ε) is equivalent to a correction of order O(u*).")
print("\nTherefore, the critical exponent ν acquires its initial non-vanishing contribution at the first order (linear) in the coupling constant u.")