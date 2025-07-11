# Step 1: Calculate the total value of Betty's estate.

# Asset 1: Flat A, worth 2,000,000 HKD.
# Betty (54) was younger than Alex (57). Under Hong Kong's Conveyancing and Property Ordinance,
# she is presumed to have survived him. The flat, held in joint tenancy, passes to her
# and becomes part of her estate.
flat_a = 2000000

# Asset 2: Flat B, worth 4,000,000 HKD, owned solely by Betty.
flat_b = 4000000

# Asset 3: Cash in bank accounts.
cash = 50000

# Asset 4: NYSE shares.
shares = 30000

# Asset 5: Personal items.
personal_items = 20000

# Asset 6: Jewellery collection.
jewellery = 500000

# Flat C is not included as it was sold before her death.

# Calculate the total estate value by summing all assets.
total_estate = flat_a + flat_b + cash + shares + personal_items + jewellery

# Step 2: Calculate the total value of valid legacies to be paid from the estate.

# Legacy 1 (Clause 4, Flat C to Cindy): Fails due to ademption (the asset was sold).
legacy_flat_c = 0

# Legacy 2 (Clause 5a, 500,000 HKD to friends): Fails due to uncertainty of beneficiaries.
legacy_friends = 0

# Legacy 3 (Clause 5b, 230,000 HKD to Wills Lawyers & Co): Valid.
legacy_lawyers = 230000

# Legacy 4 (Clause 6, 150,000 HKD to RSPCA): Fails as the clause was deleted (revoked).
legacy_rspca = 0

# Calculate the total value of legacies to be paid.
total_payouts = legacy_flat_c + legacy_friends + legacy_lawyers + legacy_rspca

# Step 3: Calculate the residuary estate.
# Residuary Estate = Total Estate - Total Payouts
residuary_estate = total_estate - total_payouts

# Print the final calculation.
print(f"Total Estate ({flat_a:,} + {flat_b:,} + {cash:,} + {shares:,} + {personal_items:,} + {jewellery:,}) = {total_estate:,} HKD")
print(f"Total Valid Legacies to be Paid Out = {total_payouts:,} HKD")
print(f"Residuary Estate Calculation: {total_estate:,} - {total_payouts:,} = {residuary_estate:,} HKD")
<<<6370000>>>