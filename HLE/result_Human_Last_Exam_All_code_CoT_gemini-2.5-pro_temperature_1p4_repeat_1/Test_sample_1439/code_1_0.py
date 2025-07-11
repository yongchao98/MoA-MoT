import sympy

# 1. Define symbolic variables for the coupling 'u', the expansion parameter 'epsilon' (ε),
# and the number of field components 'N'. The standard phi-4 theory corresponds to N=1.
u, epsilon, N = sympy.symbols('u epsilon N')

# 2. Define the one-loop Renormalization Group (RG) functions for the O(N) model in d=4-ε dimensions.
# The beta function governs the flow of the coupling constant 'u'.
beta_u = -epsilon * u + (sympy.S(N + 8) / 6) * u**2

# The anomalous dimension gamma_phi2 is related to the exponent nu.
gamma_phi2 = (sympy.S(N + 2) / 6) * u

print("--- Step 1: Defining the one-loop RG functions ---")
print("Beta function, β(u):")
sympy.pprint(beta_u)
print("\nAnomalous dimension, γ_φ²(u):")
sympy.pprint(gamma_phi2)

# 3. Find the non-trivial Wilson-Fisher fixed point, u_star, by solving beta(u) = 0.
# The solutions are u=0 (trivial Gaussian fixed point) and the non-trivial one.
u_star_solutions = sympy.solve(beta_u, u)
u_star = u_star_solutions[1]  # Select the non-trivial solution (u != 0).

print("\n--- Step 2: Finding the non-trivial fixed point u* ---")
print("Solving β(u) = 0 gives the fixed point u*:")
sympy.pprint(sympy.Eq(sympy.Symbol('u*'), u_star))
print("\nNote that the fixed point value u* is of the first order in ε.")

# 4. Define the critical exponent nu. Its inverse is given by 2 - gamma_phi2.
# First, let's look at nu as a function of u before evaluating at the fixed point.
nu_of_u = 1 / (2 - gamma_phi2)

# Expanding nu(u) in powers of u shows the structure of contributions.
nu_u_series = nu_of_u.series(u, 0, 2).removeO()

print("\n--- Step 3: Analyzing the structure of ν(u) ---")
print("The exponent ν is given by ν = 1 / (2 - γ_φ²(u)).")
print("Expanding ν as a series in the coupling constant u:")
sympy.pprint(sympy.Eq(sympy.Symbol('ν(u)'), nu_u_series))
print("\nThis shows that the first correction to the mean-field value (1/2) is of the first order in u, i.e., O(u).")

# 5. Substitute the fixed point u_star into the expression for nu to get its value at criticality.
nu_at_fixed_point = nu_of_u.subs(u, u_star)

print("\n--- Step 4: Calculating ν at the fixed point u* ---")
print("Substituting u = u* into ν(u) gives:")
sympy.pprint(sympy.Eq(sympy.Symbol('ν'), nu_at_fixed_point))

# 6. Perform a series expansion in epsilon to get the well-known result for nu.
nu_expansion = nu_at_fixed_point.series(epsilon, 0, 2).removeO()
print("\nExpanding this expression for small ε gives:")
sympy.pprint(sympy.Eq(sympy.Symbol('ν'), nu_expansion))

# 7. Print the final conclusion and the equation with numbers.
print("\n--- Conclusion ---")
print("The critical exponent ν acquires its initial non-vanishing contribution (the correction to the mean-field value of 1/2) from the term that is first-order in the coupling constant u.")
print("This corresponds to the one-loop contribution in perturbation theory.")
print("\nThe final equation for the exponent ν to first order in ε is:")
# The following line explicitly prints each number/term in the equation.
nu_coeff_numerator = N + 2
nu_coeff_denominator = 4 * (N + 8)
print(f"ν = 1/2 + ( ({nu_coeff_numerator}) / ({nu_coeff_denominator}) ) * ε + O(ε²)")
print(f"The order of the first contribution in the coupling u is therefore 1.")
