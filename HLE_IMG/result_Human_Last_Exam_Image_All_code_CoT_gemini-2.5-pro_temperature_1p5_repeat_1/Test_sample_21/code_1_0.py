# Let's model the change in the economic variables.
# We'll represent the increase in government spending with an arbitrary value, say 100 units.
delta_G = 100

print(f"1. A large fiscal expansion is implemented.")
print(f"   - Change in Government Spending (ΔG) = +{delta_G}")
print("\n")

# 2. The expansion is debt-financed, meaning it's funded by borrowing, not new taxes.
# Public Savings (Sg) = Taxes (T) - Government Spending (G).
# Since T is constant and G increases by 100, public savings decrease by 100.
delta_Sg = -delta_G
print(f"2. The expansion is debt-financed, so public savings decrease.")
print(f"   - Change in Public Savings (ΔSg) = {delta_Sg}")
print("\n")

# 3. Increased government spending is completely offset by increased private saving.
# This is the key assumption (related to Ricardian Equivalence).
delta_Sp = delta_G
print(f"3. Private savings increase to offset the government deficit.")
print(f"   - Change in Private Savings (ΔSp) = +{delta_Sp}")
print("\n")

# 4. Calculate the change in National Savings (S).
# National Savings (S) = Private Savings (Sp) + Public Savings (Sg).
# Change in National Savings (ΔS) = ΔSp + ΔSg
delta_S = delta_Sp + delta_Sg
print(f"4. Calculate the total change in National Savings (ΔS).")
print(f"   ΔS = ΔSp + ΔSg")
print(f"   ΔS = {delta_Sp} + ({delta_Sg})")
print(f"   ΔS = {delta_S}")
print("\n")

# 5. Determine the impact on the Current Account (CA).
# The Current Account is the difference between National Savings and Investment (I).
# CA = S - I
# We assume no change in investment, so ΔI = 0.
delta_I = 0
# The change in the Current Account (ΔCA) = ΔS - ΔI.
delta_CA = delta_S - delta_I
print(f"5. Calculate the final impact on the Current Account (ΔCA).")
print(f"   ΔCA = ΔS - ΔI")
print(f"   ΔCA = {delta_S} - {delta_I}")
print(f"   ΔCA = {delta_CA}")
print("\n")

print("Conclusion: The decrease in public saving is perfectly matched by an increase in private saving.")
print("This results in no change to national saving, and therefore, no change to the current account balance.")
