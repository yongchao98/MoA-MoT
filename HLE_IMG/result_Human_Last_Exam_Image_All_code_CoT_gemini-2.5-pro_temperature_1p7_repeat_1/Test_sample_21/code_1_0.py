import sys

# Define a hypothetical value for the increase in government spending (fiscal expansion).
# Any positive number will demonstrate the principle.
delta_G = 100.0

print("This program demonstrates the theoretical impact of a specific fiscal policy on the current account.")
print("-" * 70)
print("The core macroeconomic identity is: Current Account (CA) = National Saving (S) - Investment (I)")
print("We need to determine the change in each component (ΔCA, ΔS, ΔI) to find the answer.")
print("-" * 70)

# Step 1: Model the change in Government Saving (Sg)
# A debt-financed fiscal expansion means G increases and T (taxes) are constant.
# Sg = T - G. Therefore, ΔSg = ΔT - ΔG. Since ΔT=0, ΔSg = -ΔG.
delta_Sg = -delta_G
print(f"1. A debt-financed fiscal expansion of {delta_G} units occurs.")
print(f"   - Change in Government Spending (ΔG) = {delta_G}")
print(f"   - The change in government saving (ΔSg) is the negative of this amount.")
print(f"   - ΔSg = -{delta_G} = {delta_Sg}")
print("-" * 70)

# Step 2: Model the change in Private Saving (Sp)
# The problem states this is completely offset by increased private saving. This is Ricardian Equivalence.
# ΔSp = ΔG
delta_Sp = delta_G
print(f"2. Increased private saving completely offsets the fiscal expansion.")
print(f"   - This means the change in private saving (ΔSp) equals the change in government spending.")
print(f"   - ΔSp = {delta_G}")
print("-" * 70)

# Step 3: Calculate the change in National Saving (S)
# ΔS = ΔSp + ΔSg
delta_S = delta_Sp + delta_Sg
print("3. Calculate the total change in National Saving (ΔS).")
print("   - Equation: ΔS = ΔSp + ΔSg")
print(f"   - Plugging in the values: ΔS = {delta_Sp} + ({delta_Sg})")
print(f"   - Result: ΔS = {delta_S}")
print("   - National Saving does not change.")
print("-" * 70)


# Step 4: Determine the change in Investment (I)
# Since National Saving (the supply of loanable funds) is unchanged, investment is also unchanged.
delta_I = 0.0
print("4. Determine the change in Investment (ΔI).")
print(f"   - Since National Saving is unchanged (ΔS = {delta_S}), the supply of loanable funds does not shift.")
print(f"   - This implies the interest rate and the level of investment are also unchanged.")
print(f"   - ΔI = {delta_I}")
print("-" * 70)

# Step 5: Calculate the final change in the Current Account (CA)
# ΔCA = ΔS - ΔI
delta_CA = delta_S - delta_I
print("5. Calculate the final change in the Current Account (ΔCA).")
print("   - Equation: ΔCA = ΔS - ΔI")
print(f"   - Plugging in the values: ΔCA = {delta_S} - {delta_I}")
print(f"   - Result: ΔCA = {delta_CA}")
print("-" * 70)

# Final Conclusion
print("\nConclusion: The fiscal expansion has no impact on the country's current account balance.")
print("This is because the decrease in government saving is perfectly matched by an increase in private saving, leaving national saving and investment unchanged.")
# The actual final answer will be returned in the requested format below.
# The following line is a placeholder to be captured for the final answer.
# print("<<<The current account balance will not change.>>>")
