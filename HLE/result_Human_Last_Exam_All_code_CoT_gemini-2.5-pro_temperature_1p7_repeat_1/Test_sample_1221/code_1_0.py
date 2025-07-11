#
# This script calculates the value of Betty's residuary estate based on the provided will and asset information.
#

# --- Step 1: Calculate the Total Gross Estate ---

# Flat A (Orchid Apartment Complex): 2,000,000 HKD
# Under Hong Kong's commorientes rule, Betty (younger) is presumed to have survived Alex (elder).
# The property, owned in joint tenancy, passes to Betty by survivorship and becomes part of her estate.
flat_a_value = 2000000

# Flat B (Peach Apartment Complex): 4,000,000 HKD
# Owned solely by Betty.
flat_b_value = 4000000

# Cash in bank accounts: 50,000 HKD
cash_value = 50000

# Shares on NYSE: 30,000 HKD
shares_value = 30000

# Personal items: 20,000 HKD
personal_items_value = 20000

# Jewellery collection: 500,000 HKD
jewellery_value = 500000

# Summing up all assets to find the total gross estate
total_gross_estate = flat_a_value + flat_b_value + cash_value + shares_value + personal_items_value + jewellery_value

# --- Step 2: Identify and Subtract Valid Deductions ---

# Gift of Flat C to Cindy (Clause 4): 0 HKD
# The gift fails by 'ademption' as the property was sold.
gift_cindy = 0

# Legacy to Friends (Clause 5a): 0 HKD
# The gift is 'void for uncertainty' as the schedule of friends was never created.
legacy_friends = 0

# Legacy to Wills Lawyers & Co (Clause 5b): 230,000 HKD
# This is a valid pecuniary legacy.
legacy_lawyers = 230000

# Legacy to RSPCA (Clause 6): 0 HKD
# This clause was deleted, thus the gift was 'revoked'.
legacy_rspca = 0

# Debts and funeral/administrative expenses are not specified and are assumed to be 0 for this calculation.
other_liabilities = 0

# Summing up all valid deductions
total_deductions = gift_cindy + legacy_friends + legacy_lawyers + legacy_rspca + other_liabilities

# --- Step 3: Calculate the Residuary Estate ---
residuary_estate = total_gross_estate - total_deductions

# --- Final Output ---
print("Calculation of Betty's Residuary Estate (in HKD)")
print("==================================================")
print("\nPart 1: Total Gross Estate Calculation")
print(f"  Flat A (by survivorship):       {flat_a_value:>12,}")
print(f"  Flat B:                         {flat_b_value:>12,}")
print(f"  Cash in Bank:                   {cash_value:>12,}")
print(f"  Shares:                         {shares_value:>12,}")
print(f"  Personal Items:                 {personal_items_value:>12,}")
print(f"  Jewellery:                      {jewellery_value:>12,}")
print("--------------------------------------------------")
print(f"  Total Gross Estate:             {total_gross_estate:>12,}")
print("\nPart 2: Total Deductions Calculation")
print(f"  Legacy to Wills Lawyers & Co:   {legacy_lawyers:>12,}")
print("  (Other gifts failed or were revoked)")
print("--------------------------------------------------")
print(f"  Total Deductions:               {total_deductions:>12,}")

print("\n==================================================")
print("Final Residuary Estate Calculation:")
print(f"{total_gross_estate:,} (Total Gross Estate) - {total_deductions:,} (Total Deductions) = {residuary_estate:,} HKD")
print("==================================================")