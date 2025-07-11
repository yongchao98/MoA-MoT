import sympy

# Define d as a symbol
d = sympy.Symbol('d')

# --- Case 1: a_1 = 2d ---
print("--- Analyzing Case 1 (a_1 = 2d) ---")
# The derived equation is 51*d^2 - d - 50 = 0
eq1 = 51*d**2 - d - 50
solutions1 = sympy.solve(eq1, d)
print("The equation is: 51*d**2 - 1*d - 50 = 0")
print(f"Potential solutions for d: {solutions1}")

# Filter for solutions where d > 1
valid_solutions1 = [sol.evalf() for sol in solutions1 if sol > 1]
print(f"Valid solutions for d > 1: {valid_solutions1}\n")


# --- Case 2: a_1 = d ---
print("--- Analyzing Case 2 (a_1 = d) ---")
# The derived equation is 50*d^2 - d - 51 = 0
eq2 = 50*d**2 - d - 51
solutions2 = sympy.solve(eq2, d)
print("The equation is: 50*d**2 - 1*d - 51 = 0")
print(f"Potential solutions for d: {solutions2}")

# Filter for solutions where d > 1
valid_solutions2 = [sol.evalf() for sol in solutions2 if sol > 1]
final_answer = None
if valid_solutions2:
    final_answer = valid_solutions2[0]
    print(f"Valid solutions for d > 1: {valid_solutions2}")
    print(f"\nThe final answer is the valid solution from Case 2.")
else:
     print("No valid solution found in Case 2 that satisfies d > 1.")


# Print the final answer
if final_answer is not None:
    print("\nFinal calculation verification:")
    # Using the equation from Case 2: 50*d - 51/d = 1
    d_val = float(final_answer)
    S99 = 4950 * d_val
    T99 = 5049 / d_val
    diff = S99 - T99
    print(f"d = {d_val}")
    print(f"S_99 = 4950 * {d_val:.2f} = {S99:.2f}")
    print(f"T_99 = 5049 / {d_val:.2f} = {T99:.2f}")
    print(f"S_99 - T_99 = {S99:.2f} - {T99:.2f} = {diff:.2f}")
    print("\nThe result matches the condition S_99 - T_99 = 99.")


<<<1.02>>>