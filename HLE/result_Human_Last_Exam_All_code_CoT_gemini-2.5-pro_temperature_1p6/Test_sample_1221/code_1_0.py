# Step 1: Define and calculate the total value of Betty's estate
# Flat A is included due to the 'commorientes rule' (right of survivorship)
# where Betty (younger) is presumed to have outlived Alex (elder).
flat_a = 2000000
# Flat B was solely owned by Betty.
flat_b = 4000000
# Other assets owned by Betty.
bank_accounts = 50000
shares = 30000
personal_items = 20000
jewellery = 500000

# Summing up all assets to find the total estate value
total_estate_value = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery

# Step 2: Define and calculate the value of valid gifts from the will
# Gift of Flat C fails due to ademption (property was sold).
gift_flat_c = 0
# Gift to friends fails for uncertainty (beneficiary list was never created).
gift_to_friends = 0
# Gift to the law firm is a valid pecuniary legacy.
gift_to_lawyers = 230000
# Gift to RSPCA fails as the clause was revoked (deleted).
gift_to_rspca = 0

# Summing up all valid gifts (legacies) to be paid out
total_valid_gifts = gift_flat_c + gift_to_friends + gift_to_lawyers + gift_to_rspca

# Step 3: Calculate the residuary estate
residuary_estate = total_estate_value - total_valid_gifts

# Print the final calculation showing each component
print("Calculation of Betty's Residuary Estate in HKD")
print("-" * 50)
print(f"Total Estate Value: {flat_a} (Flat A) + {flat_b} (Flat B) + {bank_accounts} (Bank) + {shares} (Shares) + {personal_items} (Items) + {jewellery} (Jewellery) = {total_estate_value}")
print(f"Total Valid Gifts: {gift_flat_c} (Flat C) + {gift_to_friends} (Friends) + {gift_to_lawyers} (Lawyers) + {gift_to_rspca} (RSPCA) = {total_valid_gifts}")
print("-" * 50)
print(f"Residuary Estate = Total Estate Value ({total_estate_value}) - Total Valid Gifts ({total_valid_gifts})")
print(f"Final Equation: {total_estate_value} - {total_valid_gifts} = {residuary_estate}")
print(f"\nThe total value of Betty's residuary estate is: {residuary_estate} HKD.")
