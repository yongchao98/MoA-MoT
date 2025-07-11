# --- User-provided variables ---
# q: Total number of outstanding shares
# E: Total market value of equity before the policy change
# d: Total dividends planned for year 1 under the original policy
# g: Annual dividend growth rate (as a decimal)

q = 10000000.0
E = 200000000.0
d = 8000000.0
g = 0.04

# --- Calculation ---
# The formula derived is p1 = (E * (1 + g) - d) / q
p1 = (E * (1 + g) - d) / q

# --- Output ---
print("The formula for the per-share ex-dividend price in year 1 (p1) is:")
print("p1 = (E * (1 + g) - d) / q")
print("\nSubstituting the given values into the formula:")
print(f"p1 = ({int(E)} * (1 + {g}) - {int(d)}) / {int(q)}")
print(f"p1 = ({E * (1 + g)} - {int(d)}) / {int(q)}")
print(f"p1 = {E * (1 + g) - d} / {int(q)}")
print(f"\nThe final per-share ex-dividend price (p1) is: ${p1:.2f}")
