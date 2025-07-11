# Step 1: Define the value of assets forming Betty's estate.

# Flat A is included fully due to the commorientes rule. Betty (54) is presumed to survive Alex (57).
flat_a = 2000000
# Flat B was solely owned by Betty.
flat_b = 4000000
# Cash, shares, and personal items are part of the estate.
cash_in_bank = 50000
nyse_shares = 30000
personal_items = 20000
# Jewellery collection is part of the estate.
jewellery_collection = 500000

# The gift of Flat C fails due to ademption (it was sold). So, its value is not added.

# Calculate the total gross value of the estate.
total_assets = flat_a + flat_b + cash_in_bank + nyse_shares + personal_items + jewellery_collection

print("Calculating Total Gross Estate:")
print(f"  Flat A (via survivorship): {flat_a}")
print(f"  Flat B (sole ownership): {flat_b}")
print(f"  Cash in Bank: {cash_in_bank}")
print(f"  NYSE Shares: {nyse_shares}")
print(f"  Personal Items: {personal_items}")
print(f"  Jewellery Collection: {jewellery_collection}")
print(f"Total Gross Estate Value = {flat_a} + {flat_b} + {cash_in_bank} + {nyse_shares} + {personal_items} + {jewellery_collection} = {total_assets} HKD")
print("-" * 30)

# Step 2: Define the value of valid legacies to be deducted.

# Gift to Cindy (Flat C) failed (adeemed).
legacy_cindy = 0
# Gift to friends failed (uncertainty of objects).
legacy_friends = 0
# Gift to Wills Lawyers & Co is valid.
legacy_lawyers = 230000
# Gift to RSPCA was revoked (deleted clause).
legacy_rspca = 0

# Calculate the total deductions for legacies.
total_deductions = legacy_cindy + legacy_friends + legacy_lawyers + legacy_rspca

print("Calculating Total Deductions (Valid Legacies):")
print(f"  Legacy to Cindy: {legacy_cindy} (Failed by ademption)")
print(f"  Legacy to Friends: {legacy_friends} (Failed by uncertainty)")
print(f"  Legacy to Wills Lawyers & Co: {legacy_lawyers} (Valid)")
print(f"  Legacy to RSPCA: {legacy_rspca} (Failed by revocation)")
print(f"Total Deductions = {legacy_cindy} + {legacy_friends} + {legacy_lawyers} + {legacy_rspca} = {total_deductions} HKD")
print("-" * 30)

# Step 3: Calculate the residuary estate.
# Residuary Estate = Total Gross Estate - Total Deductions
# Note: This is before considering debts and administrative/funeral expenses.

residuary_estate = total_assets - total_deductions

print("Final Calculation of the Residuary Estate:")
print(f"Residuary Estate = Total Assets - Total Deductions")
print(f"Residuary Estate = {total_assets} - {total_deductions}")
print(f"Total Value of Residuary Estate = {residuary_estate} HKD")
