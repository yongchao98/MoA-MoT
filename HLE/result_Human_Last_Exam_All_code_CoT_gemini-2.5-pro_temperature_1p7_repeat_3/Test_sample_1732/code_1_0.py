import sympy

# Define d as a symbol to solve for it
d = sympy.symbols('d')

print("Step 1: Determine the relationship between a1 and d.")
print("Since {b_n} is an arithmetic progression, we have 2*b_2 = b_1 + b_3.")
print("Substituting b_n = n*(n+1) / (a1 + (n-1)*d) leads to the equation:")
print("a1^2 - 3*a1*d + 2*d^2 = 0")
print("This factors to (a1 - d)*(a1 - 2*d) = 0, giving two possible cases.")
print("Case 1: a1 = d")
print("Case 2: a1 = 2*d\n")

# --- Case 1: a1 = d ---
print("--- Analyzing Case 1: a1 = d ---")
print("If a1 = d, then a_n = d + (n-1)*d = n*d.")
print("And b_n = n*(n+1) / (n*d) = (n+1)/d.")
N = 99

# S_99 is the sum of a_n = n*d from n=1 to 99
# S_99 = d * (1 + 2 + ... + 99) = d * (99 * 100 / 2)
S99_val_1 = 4950 * d

# T_99 is the sum of b_n = (n+1)/d from n=1 to 99
# T_99 = (1/d) * (2 + 3 + ... + 100) = (1/d) * ((100 * 101 / 2) - 1)
T99_val_1 = 5049 / d

print(f"S_99 = {S99_val_1}")
print(f"T_99 = {T99_val_1}")

# The final equation using S_99 - T_99 = 99
# 4950*d - 5049/d = 99
# Dividing by 99 gives 50*d - 51/d = 1
# This results in the quadratic equation 50*d^2 - d - 51 = 0
s1_coeff = 4950
t1_coeff = 5049
const = 99
print(f"\nThe equation S_99 - T_99 = 99 becomes: {s1_coeff}*d - {t1_coeff}/d = {const}")

quad_c1 = s1_coeff // const
quad_c2 = const // const
quad_c3 = t1_coeff // const
print(f"Rearranging this into a quadratic equation gives: {quad_c1}*d**2 - {quad_c2}*d - {quad_c3} = 0")

solutions_1 = sympy.solve(50*d**2 - d - 51, d)
print(f"The solutions for d are {solutions_1[0]} and {solutions_1[1]}.")
valid_d = None
for sol in solutions_1:
    if sol > 1:
        valid_d = sol
        print(f"The solution d = {sol} (or {float(sol)}) satisfies the condition d > 1.\n")

# --- Case 2: a1 = 2d ---
print("--- Analyzing Case 2: a1 = 2d ---")
print("If a1 = 2d, then a_n = 2d + (n-1)*d = (n+1)*d.")
print("And b_n = n*(n+1) / ((n+1)*d) = n/d.")

# S_99 is the sum of a_n = (n+1)*d from n=1 to 99
# S_99 = d * (2 + 3 + ... + 100)
S99_val_2 = 5049 * d

# T_99 is the sum of b_n = n/d from n=1 to 99
# T_99 = (1/d) * (1 + 2 + ... + 99)
T99_val_2 = 4950 / d

print(f"S_99 = {S99_val_2}")
print(f"T_99 = {T99_val_2}")

s2_coeff = 5049
t2_coeff = 4950
print(f"\nThe equation S_99 - T_99 = 99 becomes: {s2_coeff}*d - {t2_coeff}/d = {const}")

quad_c1_2 = s2_coeff // const
quad_c2_2 = const // const
quad_c3_2 = t2_coeff // const
print(f"Rearranging this into a quadratic equation gives: {quad_c1_2}*d**2 - {quad_c2_2}*d - {quad_c3_2} = 0")

solutions_2 = sympy.solve(51*d**2 - d - 50, d)
print(f"The solutions for d are {solutions_2[0]} and {solutions_2[1]}.")
print("Neither of these solutions is greater than 1.\n")

print("--- Conclusion ---")
print(f"The only valid value for d is the one found in Case 1.")
print(f"The final answer for d is {float(valid_d)}.")
