import sympy

# Define the variable for integration
x = sympy.Symbol('x')

# The integrand is Phi(T1(sqrt(2)*x) + T2(x))
# From our analysis, T1(sqrt(2)*x) + T2(x) simplifies to x - 1.
Z_x = x - 1

# The heat flux function is Phi(z) = 0.5*z + 1.
# So the integrand is 0.5*(x - 1) + 1.
integrand = 0.5 * Z_x + 1

# Calculate the definite integral from 0 to 1.
total_energy = sympy.integrate(integrand, (x, 0, 1))

# We can also express the calculation step-by-step for clarity.
term1_coeff = 0.5
term2_const = 0.5
integral_term1 = sympy.integrate(term1_coeff * x, (x, 0, 1))
integral_term2 = sympy.integrate(term2_const, (x, 0, 1))

print("The problem asks for the minimum value of the total heat energy, given by the integral of Phi(T1 + T2).")
print("Step 1: The heat flux function Phi(z) is determined to be 0.5*z + 1.")
print("Step 2: The sum of the temperature profiles T1(sqrt(2)*x) + T2(x) simplifies to x - 1.")
print("Step 3: The expression for the total energy E is the integral from 0 to 1 of (0.5 * (x - 1) + 1), which simplifies to the integral of (0.5*x + 0.5).")
print(f"Step 4: We calculate the integral of the first term, 0.5*x, from 0 to 1, which is {integral_term1}.")
print(f"Step 5: We calculate the integral of the second term, 0.5, from 0 to 1, which is {integral_term2}.")
print(f"Step 6: The total energy is the sum of these values.")
print(f"E_total = {integral_term1} + {integral_term2} = {total_energy}")
print("Final Equation:")
print(f"{integral_term1} + {integral_term2} = {total_energy}")
