def print_calculation_step(step_number, description, equation_str):
    """Helper function to format and print each step of the calculation."""
    print(f"Step {step_number}: {description}")
    print(f"    {equation_str}\n")

print("--- Titan Computer Calculation for Firework Force F ---\n")
print("The governing equation for the force is F = 2 * m * g * sqrt(2).")
print("With mass m = 0.15 * pi kg, this becomes F = 0.3 * pi * g * sqrt(2).\n")
print("We will now compute this using 5-bit fractional arithmetic.\n")

# Initial formula with the constant as a fraction
print_calculation_step(1, "Represent the formula with fractions.", "F = (3/10) * g * pi * sqrt(2)")

# Step 2: Substitute fractional approximations for the physical constants.
# We choose values that will simplify well in later steps.
# g = 9.8 ≈ 10/1
# pi = 3.14159... ≈ 22/7
# sqrt(2) = 1.414... ≈ 7/5
g_num, g_den = 10, 1
pi_num, pi_den = 22, 7
sqrt2_num, sqrt2_den = 7, 5
print_calculation_step(2, "Approximate constants g, pi, and sqrt(2).", f"F = (3/10) * ({g_num}/{g_den}) * ({pi_num}/{pi_den}) * ({sqrt2_num}/{sqrt2_den})")

# Step 3: Reorder and perform the first multiplications.
# We can multiply (3/10) by (10/1) and (22/7) by (7/5) because they allow for cancellation.
# (3/10) * (10/1) = 3/1
# (22/7) * (7/5) = 22/5
term1_num, term1_den = 3, 1
term2_num, term2_den = 22, 5
print_calculation_step(3, "Simplify pairs of fractions by cancellation.", f"F = ({term1_num}/{term1_den}) * ({term2_num}/{term2_den})")

# Step 4: The next multiplication (3 * 22) would be 66, which exceeds the 31 limit.
# We must use the expansion strategy.
print("Step 4: The operation (3/1) * (22/5) is invalid as 3 * 22 = 66, which exceeds the 5-bit limit (31).")
print("    We must expand the fraction 22/5.\n")

# Step 5: Expand 22/5 into a sum of valid fractions.
# 22/5 = 4 + 2/5 = 4/1 + 2/5
expanded_term1_num, expanded_term1_den = 4, 1
expanded_term2_num, expanded_term2_den = 2, 5
print_calculation_step(5, "Expand 22/5 into (4/1 + 2/5).", f"F = ({term1_num}/{term1_den}) * (({expanded_term1_num}/{expanded_term1_den}) + ({expanded_term2_num}/{expanded_term2_den}))")

# Step 6: Distribute the multiplication.
# F = (3/1 * 4/1) + (3/1 * 2/5) = 12/1 + 6/5
sum_term1_num, sum_term1_den = 12, 1
sum_term2_num, sum_term2_den = 6, 5
print_calculation_step(6, "Distribute the multiplication.", f"F = ({sum_term1_num}/{sum_term1_den}) + ({sum_term2_num}/{sum_term2_den})")

# Step 7: The addition requires a common denominator, which would be (60/5) + (6/5).
# The numerator 60 exceeds the 31 limit. We must approximate the smaller term.
print("Step 7: The operation (12/1) + (6/5) is invalid as the common denominator form (60/5) + (6/5) has a numerator of 60, exceeding 31.")
print("    To proceed, we approximate the smaller term, 6/5 (value 1.2), with the closest valid fraction, 1/1.\n")

# Step 8: Replace 6/5 with its approximation 1/1 and perform the final addition.
approx_term_num, approx_term_den = 1, 1
final_num = sum_term1_num * approx_term_den + approx_term_num * sum_term1_den
final_den = sum_term1_den * approx_term_den
print_calculation_step(8, "Approximate 6/5 with 1/1 and add.", f"F = ({sum_term1_num}/{sum_term1_den}) + ({approx_term_num}/{approx_term_den}) = {final_num}/{final_den}")

print(f"--- Final Result ---\nThe calculated force F is {final_num}/{final_den}, which is {float(final_num/final_den)} N.")
target_F = 0.3 * 9.8 * 3.14159265 * 1.41421356
error = abs(final_num - target_F)
print(f"The precise target value is ~{target_F:.4f} N.")
print(f"The absolute error is |{final_num} - {target_F:.4f}| = {error:.4f} N.")