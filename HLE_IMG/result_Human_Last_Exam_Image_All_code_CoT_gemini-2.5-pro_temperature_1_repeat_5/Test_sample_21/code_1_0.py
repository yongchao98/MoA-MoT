# Assume a hypothetical fiscal expansion amount for demonstration.
# Let the increase in government spending (G) be 100 units.
delta_G = 100

# According to the problem statement:
# 1. The fiscal expansion is debt-financed, so public saving falls by the amount of spending.
delta_S_public = -delta_G

# 2. The increased government spending is completely offset by increased private saving.
delta_S_private = delta_G

# The total change in national saving (S) is the sum of the changes in private and public saving.
delta_S = delta_S_private + delta_S_public

# Since national saving does not change, the interest rate is unaffected, and thus investment (I) does not change.
delta_I = 0

# The change in the Current Account (CA) is the change in National Saving minus the change in Investment.
delta_CA = delta_S - delta_I

print("The relationship is: ΔCA = ΔS_private + ΔS_public - ΔI")
print("Substituting the values based on the scenario:")
print(f"ΔCA = {delta_S_private} + ({delta_S_public}) - {delta_I}")
print(f"ΔCA = {delta_CA}")
print("\nConclusion: The fiscal expansion has no impact on the country's current account balance.")