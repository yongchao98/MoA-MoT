# Titan Superconducting Computer Simulation

# This code demonstrates the final crucial calculation step for the force F,
# showing how approximation is used to overcome the 5-bit integer overflow.

# Initial values after several steps of simplification
g_term_num = 10
g_term_den = 1
other_term_num = 9
other_term_den = 7

# The multiplication (g_term_num * other_term_num) = 10 * 9 = 90, which exceeds 31.
# We must approximate one of the terms.
# We approximate g = 10/1 with g' = 21/2 (10.5) to allow for simplification.
g_approx_num = 21
g_approx_den = 2

print("Problem: Cannot compute (10/1) * (9/7) directly due to 10 * 9 = 90 > 31.")
print("Solution: Approximate g=10/1 with g'=21/2.")
print("New Calculation: (21/2) * (9/7)")

# Simplify before multiplying: 21 is divisible by 7
common_divisor = 7
simplified_g_num = g_approx_num // common_divisor
simplified_other_den = other_term_den // common_divisor
print(f"Simplify by dividing {g_approx_num} and {other_term_den} by {common_divisor}.")
print(f"Expression becomes: ({simplified_g_num}/{g_approx_den}) * ({other_term_num}/{simplified_other_den})")

# Final multiplication
final_num = simplified_g_num * other_term_num
final_den = g_approx_den * simplified_other_den
print(f"Final multiplication: ({simplified_g_num} * {other_term_num}) / ({g_approx_den} * {simplified_other_den}) = {final_num}/{final_den}")

# Calculate true force for error analysis
import math
F_true = (9.8 * 20 * (0.15 * math.pi)) / (10 * math.cos(math.radians(45)))
F_calc = final_num / final_den
error = abs(F_true - F_calc)

print(f"\nCalculated Force F = {F_calc}")
print(f"True Force F_true = {F_true:.3f}")
print(f"Absolute Error e = |{F_true:.3f} - {F_calc}| = {error:.3f}")

# Final Answer format
final_answer = f"Y[{error:.3f}]"
print("\nFinal Answer:")
print(f"<<<{final_answer}>>>")