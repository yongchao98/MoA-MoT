# This script calculates the value of Betty's residuary estate based on Hong Kong law.

# 1. Define the value of assets in Betty's estate.
# Flat A is included due to the commorientes rule (younger survives elder in simultaneous deaths).
flat_a_value = 2000000
flat_b_value = 4000000
bank_accounts = 50000
shares = 30000
personal_items = 20000
jewellery = 500000

# Calculate the total value of the estate.
total_estate = flat_a_value + flat_b_value + bank_accounts + shares + personal_items + jewellery

# 2. Define the value of valid legacies to be paid from the estate.
# Gift of Flat C to Cindy failed (adeemed).
legacy_cindy_flat_c = 0
# Gift to friends failed (uncertainty of objects).
legacy_friends_schedule = 0
# Gift to Wills Lawyers & Co is valid.
legacy_wills_lawyers = 230000
# Gift to RSPCA failed (revoked).
legacy_rspca = 0

# Calculate the total value of legacies to be paid out.
total_legacies = legacy_cindy_flat_c + legacy_friends_schedule + legacy_wills_lawyers + legacy_rspca

# 3. Calculate the residuary estate.
residuary_estate = total_estate - total_legacies

# 4. Print the final calculation and result.
print("Calculation of Betty's Residuary Estate:")
print(f"Total Estate Value = {flat_a_value} (Flat A) + {flat_b_value} (Flat B) + {bank_accounts} (Bank) + {shares} (Shares) + {personal_items} (Items) + {jewellery} (Jewellery) = {total_estate} HKD")
print(f"Total Valid Legacies = {legacy_wills_lawyers} (to Wills Lawyers & Co) = {total_legacies} HKD")
print("\nFinal Equation:")
print(f"{total_estate} (Total Estate) - {total_legacies} (Valid Legacies) = {residuary_estate} (Residuary Estate)")
print(f"\nThe total value of Betty's residuary estate is {residuary_estate} HKD.")
