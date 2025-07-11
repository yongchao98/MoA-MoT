# Goal: To determine the theoretical impact of a specific fiscal expansion on the current account.

# 1. Define the fundamental identity.
print("Step 1: The current account (CA) is defined as National Saving (S) minus Investment (I).")
print("CA = S - I\n")

# 2. Define the components of National Saving.
print("Step 2: National Saving (S) is the sum of Private Saving (Sp) and Public Saving (Sg).")
print("S = Sp + Sg\n")

# 3. Model the policy shock. Let's use a hypothetical value for the fiscal expansion.
# Assume the government increases spending by $100 billion.
fiscal_expansion_G = 100

# The fiscal expansion is debt-financed, meaning it reduces public saving by the full amount.
change_in_public_saving_Sg = -fiscal_expansion_G

# The problem states that private saving increases to completely offset the government spending.
change_in_private_saving_Sp = fiscal_expansion_G

print(f"Step 3: Analyze the policy's impact on saving components.")
print(f"The fiscal expansion means government spending increases by {fiscal_expansion_G}.")
print(f"This reduces Public Saving by {fiscal_expansion_G}, so the change in Sg is {change_in_public_saving_Sg}.")
print(f"Private Saving increases to offset this, so the change in Sp is {change_in_private_saving_Sp}.\n")


# 4. Calculate the total change in National Saving.
change_in_national_saving_S = change_in_private_saving_Sp + change_in_public_saving_Sg

print("Step 4: Calculate the total change in National Saving (S).")
print(f"Change in S = (Change in Sp) + (Change in Sg)")
print(f"Change in S = {change_in_private_saving_Sp} + ({change_in_public_saving_Sg}) = {change_in_national_saving_S}\n")
print("Conclusion for Step 4: There is no change in total National Saving.\n")

# 5. Assess the impact on Investment (I).
# In this theoretical framework (Ricardian Equivalence), since national saving is unchanged,
# there is no pressure on interest rates, so investment is also unchanged.
change_in_investment_I = 0

print(f"Step 5: Assume that with no change in national saving, Investment (I) also does not change.")
print(f"Change in I = {change_in_investment_I}\n")


# 6. Calculate the final change in the Current Account.
change_in_current_account_CA = change_in_national_saving_S - change_in_investment_I

print("Step 6: Calculate the final impact on the Current Account (CA).")
print("Change in CA = (Change in S) - (Change in I)")
# This is the final equation where each number must be outputted.
print(f"Change in CA = {change_in_national_saving_S} - {change_in_investment_I} = {change_in_current_account_CA}\n")

# Final Answer
print("Final Answer: The fiscal expansion has no impact on the country's current account balance.")