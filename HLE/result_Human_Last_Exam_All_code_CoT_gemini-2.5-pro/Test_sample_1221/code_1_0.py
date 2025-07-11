# Step 1: Define the value of assets in Betty's estate
# Flat A passed to Betty by survivorship as she was younger than Alex.
flat_a_value = 2000000
# Flat B was solely owned by Betty.
flat_b_value = 4000000
# Other personal assets.
bank_account_cash = 50000
shares_value = 30000
personal_items_value = 20000
jewellery_value = 500000

# Calculate the total value of the estate
total_estate_assets = flat_a_value + flat_b_value + bank_account_cash + shares_value + personal_items_value + jewellery_value

# Step 2: Define the value of valid legacies to be paid out
# Gift of Flat C to Cindy fails due to ademption.
legacy_flat_c = 0
# Gift to friends fails due to uncertainty of objects.
legacy_friends = 0
# Gift to RSPCA was deleted from the will.
legacy_rspca = 0
# Gift to the law firm is a valid legacy.
legacy_lawyers = 230000

# Calculate the total value of legacies to be deducted from the estate
total_legacies_paid_out = legacy_flat_c + legacy_friends + legacy_rspca + legacy_lawyers

# Step 3: Calculate the residuary estate
residuary_estate = total_estate_assets - total_legacies_paid_out

# Print the final calculation and result
print("Calculation of Betty's Residuary Estate:")
print(f"Total Assets: {flat_a_value} (Flat A) + {flat_b_value} (Flat B) + {bank_account_cash} (Cash) + {shares_value} (Shares) + {personal_items_value} (Items) + {jewellery_value} (Jewellery) = {total_estate_assets} HKD")
print(f"Total Valid Legacies: {total_legacies_paid_out} HKD")
print(f"Residuary Estate = Total Assets - Total Valid Legacies")
print(f"Residuary Estate = {total_estate_assets} - {total_legacies_paid_out} = {residuary_estate} HKD")
