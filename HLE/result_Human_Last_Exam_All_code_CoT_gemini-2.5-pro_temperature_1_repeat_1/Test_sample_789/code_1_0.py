import math

# Helper function to show a calculation step
def show_calc(op, f1_n, f1_d, f2_n, f2_d, res_n, res_d, comment=""):
    print(f"Step: ({f1_n}/{f1_d}) {op} ({f2_n}/{f2_d}) = {res_n}/{res_d} {comment}")

print("### Goal: Calculate Mass = density * (4/3) * pi * r^3 ###")
print("Rule: All numerators and denominators in any calculation step must be <= 31.\n")

# Step 1: Define initial values as 5-bit fractions
print("--- Initial Values ---")
density_n, density_d = 9, 10
four_thirds_n, four_thirds_d = 4, 3
pi_n, pi_d = 22, 7  # Most accurate 5-bit approximation
r_n, r_d = 1, 2
print(f"density = {density_n}/{density_d}, 4/3 = {four_thirds_n}/{four_thirds_d}, pi = {pi_n}/{pi_d}, r = {r_n}/{r_d}\n")


# Step 2: Calculate r^3
print("--- Calculate r^3 ---")
r_squared_n = r_n * r_n
r_squared_d = r_d * r_d
show_calc("*", r_n, r_d, r_n, r_d, r_squared_n, r_squared_d, "(r^2)")

r_cubed_n = r_squared_n * r_n
r_cubed_d = r_squared_d * r_d
show_calc("*", r_squared_n, r_squared_d, r_n, r_d, r_cubed_n, r_cubed_d, "(r^3)\n")


# Step 3: Group calculations to manage complexity and stay within 5-bit limits
print("--- Intermediate Calculations ---")

# Part A: density * (4/3)
# Direct multiplication (9*4)/(10*3) = 36/30 fails because 36 and 30 are > 31.
# We must simplify before multiplying: (9/3 * 4/10) -> (3 * 2) / 5 = 6/5
part_A_n = 6
part_A_d = 5
print("Part A = density * (4/3)")
print(f"Calculating ({density_n}/{density_d}) * ({four_thirds_n}/{four_thirds_d}) by simplifying first to stay within limits.")
show_calc("*", f"{density_n} (simplified to 3)", f"{density_d} (simplified to 5)", f"{four_thirds_n} (simplified to 2)", f"{four_thirds_d} (simplified to 1)", part_A_n, part_A_d, "(Valid: 6<=31, 5<=31)\n")


# Part B: pi * r^3
# Direct multiplication (22*1)/(7*8) = 22/56 fails because 56 > 31.
# We must simplify first: (22/8 * 1/7) -> (11/4 * 1/7) = 11/28
part_B_n = 11
part_B_d = 28
print("Part B = pi * r^3")
print(f"Calculating ({pi_n}/{pi_d}) * ({r_cubed_n}/{r_cubed_d}) by simplifying first to stay within limits.")
show_calc("*", f"{pi_n} (simplified to 11)", f"{pi_d}", f"{r_cubed_n}", f"{r_cubed_d} (simplified to 4)", part_B_n, part_B_d, "(Valid: 11<=31, 28<=31)\n")


# Step 4: Final multiplication and necessary approximation
print("--- Final Calculation ---")
print("Mass = Part A * Part B")
print(f"Attempting: ({part_A_n}/{part_A_d}) * ({part_B_n}/{part_B_d})")
print(f"Check: Numerator = {part_A_n*part_B_n}, Denominator = {part_A_d*part_B_d}")
print("This fails because 66 > 31 and 140 > 31.\n")

# As per the rules, we must replace a fraction with a less precise one.
# We replace Part B (11/28) with an approximation that makes the calculation possible.
approx_B_n = 2
approx_B_d = 5
print(f"Approximation: Part B ({part_B_n}/{part_B_d} ≈ 0.393) is replaced by {approx_B_n}/{approx_B_d} = 0.400.")
print("This is the valid fraction that minimizes error while allowing the multiplication (6*2<=31 and 5*5<=31).\n")

# Final calculation
final_mass_n = part_A_n * approx_B_n
final_mass_d = part_A_d * approx_B_d
print("Performing the final calculation with the approximation:")
print(f"Final Equation: {part_A_n} / {part_A_d} * {approx_B_n} / {approx_B_d} = {final_mass_n} / {final_mass_d}")
