# Let's represent the increase in government spending (fiscal expansion) by an arbitrary value, for example, 100 units.
delta_G = 100

# According to the problem, the increase in private saving is equal to the increase in government spending.
delta_Sp = delta_G

# The expansion is debt-financed, meaning taxes (T) don't change.
# Public Saving (Sg) = T - G.
# The change in public saving (delta_Sg) is the change in T minus the change in G.
# delta_Sg = delta_T - delta_G = 0 - delta_G
delta_Sg = -delta_G

# We assume investment (I) does not change.
delta_I = 0

# The change in the Current Account (CA) is the change in National Savings minus the change in Investment.
# National Savings = Private Savings (Sp) + Public Savings (Sg)
# So, delta_CA = delta_Sp + delta_Sg - delta_I
delta_CA = delta_Sp + delta_Sg - delta_I

print("The fundamental equation for the change in the Current Account is:")
print("ΔCA = ΔSp + ΔSg - ΔI")
print("\nBased on the problem's assumptions:")
print(f"The increase in government spending (ΔG) is {delta_G}.")
print(f"The corresponding increase in private savings (ΔSp) is {delta_Sp}.")
print(f"The corresponding decrease in public savings (ΔSg) is {delta_Sg}.")
print(f"The change in investment (ΔI) is assumed to be {delta_I}.")
print("\nPlugging these values into the equation:")
print(f"ΔCA = {delta_Sp} + ({delta_Sg}) - {delta_I}")

final_equation = f"ΔCA = {delta_Sp} {delta_Sg} - {delta_I}"
print(final_equation)

print(f"\nΔCA = {delta_CA}")
print("\nConclusion: The fiscal expansion has no impact on the country's current account balance.")
print("This scenario is a classic example of Ricardian Equivalence.")
