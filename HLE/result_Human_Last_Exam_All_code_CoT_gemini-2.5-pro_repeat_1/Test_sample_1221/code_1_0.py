# Step 1: Define the value of Betty's assets.
# Under the commorientes rule, Betty, being younger, is presumed to have survived Alex.
# Thus, the joint tenancy property (Flat A) passes to her entirely before becoming part of her estate.
flat_a_value = 2000000
flat_b_value = 4000000
bank_cash = 50000
nyse_shares = 30000
personal_items = 20000
jewellery = 500000

# Calculate the total gross estate.
total_assets = flat_a_value + flat_b_value + bank_cash + nyse_shares + personal_items + jewellery

print("Calculating the Total Gross Estate:")
print(f"  Value of Flat A (survivorship): {flat_a_value:,} HKD")
print(f"  Value of Flat B: {flat_b_value:,} HKD")
print(f"  Cash in bank: {bank_cash:,} HKD")
print(f"  NYSE Shares: {nyse_shares:,} HKD")
print(f"  Personal Items: {personal_items:,} HKD")
print(f"  Jewellery: {jewellery:,} HKD")
print(f"Total Gross Estate = {total_assets:,} HKD\n")


# Step 2: Identify valid deductions from the estate based on the will.
# Gift of Flat C to Cindy failed (adeemed).
# Gift of 500,000 HKD to friends failed (uncertainty).
# Gift of 150,000 HKD to RSPCA failed (deleted).
# The only valid legacy to be paid is to the lawyers.
legacy_to_lawyers = 230000

print("Calculating Deductions from the Estate:")
print("  Gift of Flat C to Cindy: Failed (property was sold).")
print("  Gift to friends: Failed (beneficiaries not specified).")
print("  Gift to RSPCA: Failed (clause was deleted).")
print(f"  Valid legacy to Wills Lawyers & Co: {legacy_to_lawyers:,} HKD\n")


# Step 3: Calculate the residuary estate.
# The residuary estate is what remains after all valid legacies are paid.
# Note: This is calculated before the settlement of any unspecified debts mentioned in the problem.
residuary_estate = total_assets - legacy_to_lawyers

print("Calculating the Final Residuary Estate (before settling unspecified debts):")
print("Residuary Estate = Total Gross Estate - Valid Legacies")
print(f"Residuary Estate = {flat_a_value:,} + {flat_b_value:,} + {bank_cash:,} + {nyse_shares:,} + {personal_items:,} + {jewellery:,} - {legacy_to_lawyers:,}")
print(f"Residuary Estate = {total_assets:,} - {legacy_to_lawyers:,}")
print(f"Total Residuary Estate = {residuary_estate:,} HKD")
