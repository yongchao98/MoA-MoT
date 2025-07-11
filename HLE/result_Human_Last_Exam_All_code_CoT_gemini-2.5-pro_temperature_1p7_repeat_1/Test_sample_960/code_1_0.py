# The problem can be solved analytically using martingale theory.
# The derivation in the text leads to the equation p <= 2 / (3*C_T - 1),
# where p is the probability of reaching a 50% state, and C_T is the
# conditional expectation of the martingale at the stopping time.
# A detailed analysis shows that C_T >= 5/3.
# We substitute this value to find the upper bound for p.

C_T_lower_bound = 5/3
numerator = 2
denominator = 3 * C_T_lower_bound - 1

upper_bound = numerator / denominator

# We want to show the final equation with the numbers plugged in.
# Equation: p <= 2 / (3 * (5/3) - 1)
p_le_numerator = 2
p_le_denominator_term1 = 3
p_le_denominator_term2 = "(5/3)"
p_le_denominator_term3 = 1

# Equation: p <= 2 / (5 - 1)
p_le_step2_numerator = 2
p_le_step2_denominator_term1 = 5
p_le_step2_denominator_term2 = 1

# Equation: p <= 2 / 4
p_le_step3_numerator = 2
p_le_step3_denominator = 4

# Final result: p <= 0.5
final_result = upper_bound

print("The derivation leads to the inequality: p <= 2 / (3*C_T - 1)")
print("Using the known lower bound for C_T (the conditional expectation of the martingale M_T at the stopping time), which is C_T >= 5/3:")
print(f"p <= {p_le_numerator} / ({p_le_denominator_term1} * {p_le_denominator_term2} - {p_le_denominator_term3})")
print(f"p <= {p_le_step2_numerator} / ({p_le_step2_denominator_term1} - {p_le_step2_denominator_term2})")
print(f"p <= {p_le_step3_numerator} / {p_le_step3_denominator}")
print(f"The upper bound for the probability is {final_result}.")