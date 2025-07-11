# Step 1: Define the value of Betty's assets
flat_a = 2000000 # Passed to Betty's estate due to commorientes rule
flat_b = 4000000 # Solely owned by Betty
bank_accounts = 50000
nyse_shares = 30000
personal_items = 20000
jewellery = 500000

# Calculate the total value of the estate
total_estate = flat_a + flat_b + bank_accounts + nyse_shares + personal_items + jewellery

# Step 2: Define the value of the legacies to be paid out
# Gift to Cindy (Flat C) failed due to ademption
legacy_cindy = 0
# Gift to friends failed due to uncertainty of objects
legacy_friends = 0
# Gift to Wills Lawyers & Co is a valid legacy
legacy_lawyers = 230000
# Gift to RSPCA was revoked by deletion
legacy_rspca = 0

# Calculate the total value of valid legacies
total_legacies_paid = legacy_cindy + legacy_friends + legacy_lawyers + legacy_rspca

# Step 3: Calculate the residuary estate
residuary_estate = total_estate - total_legacies_paid

# Print the final equation with all numbers
print("Calculation of Betty's Residuary Estate:")
print(f"Total Estate ({flat_a} + {flat_b} + {bank_accounts} + {nyse_shares} + {personal_items} + {jewellery}) = {total_estate}")
print(f"Total Valid Legacies to be Paid Out = {total_legacies_paid}")
print(f"Residuary Estate = {total_estate} - {total_legacies_paid} = {residuary_estate}")
print(f"\nThe value of Betty's residuary estate is {residuary_estate} HKD.")
<<<6370000>>>