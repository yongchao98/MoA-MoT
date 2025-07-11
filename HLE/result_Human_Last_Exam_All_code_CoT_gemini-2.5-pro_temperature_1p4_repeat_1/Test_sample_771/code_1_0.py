import math

# Define the parameters for the calculation
p0 = 1.0  # Initial price
r_a = 0.10 # Annual growth rate for Basket A
g_b = 0.10 # Annual linear growth trend for Basket B
t = 30   # Time in years

# --- Calculate the final value for Basket A ---
# Basket A grows with compound interest and has no volatility.
# The formula is P(t) = P0 * (1 + r)^t
final_value_a = p0 * math.pow(1 + r_a, t)

# --- Calculate the expected final value for Basket B ---
# Basket B has a linear growth trend, modeled as an arithmetic random walk.
# The formula for the expected value is E[P(t)] = P0 * (1 + g * t)
expected_final_value_b = p0 * (1 + g_b * t)

# --- Output the results ---
print("Analysis of Expected Values over 30 years:")
print("-" * 40)

# Print the calculation for Basket A
print("Basket A (Deterministic Compound Growth):")
print(f"Final Value = {p0} * (1 + {r_a})^{t}")
print(f"Final Value = {final_value_a:.4f}")
print("")

# Print the calculation for Basket B
print("Basket B (Linear Growth Trend with Volatility):")
print(f"Expected Final Value = {p0} * (1 + {g_b} * {t})")
print(f"Expected Final Value = {expected_final_value_b:.4f}")
print("-" * 40)

# --- Compare the values and determine the conclusion ---
if final_value_a > expected_final_value_b:
    conclusion = "The strategy of choosing basket A has a higher expected value than that of choosing basket B."
    answer = "A"
elif expected_final_value_b > final_value_a:
    conclusion = "The strategy of choosing basket B has a higher expected value than that of choosing basket A."
    answer = "G"
else:
    conclusion = "The strategies have the same expected value."
    answer = "None"

print(f"Comparison: {final_value_a:.4f} (A) > {expected_final_value_b:.4f} (B)")
print(f"Conclusion: {conclusion}")

# The final answer corresponds to the letter of the correct statement.
# print(f"\nFinal Answer Code: <<< {answer} >>>")
<<<A>>>