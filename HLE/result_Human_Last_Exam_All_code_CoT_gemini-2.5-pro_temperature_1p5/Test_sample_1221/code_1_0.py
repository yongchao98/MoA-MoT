#
# Task: Calculate the value of Betty's residuary estate under Hong Kong law.
#

# --- Step 1: Define and calculate the total value of Betty's estate ---

# Asset Values (in HKD)
flat_a_orchid_apt = 2_000_000
flat_b_peach_apt = 4_000_000
bank_accounts = 50_000
nyse_shares = 30_000
personal_items = 20_000
jewellery_collection = 500_000

# Explanation of Assets:
# - Flat A: Included in full due to the right of survivorship. Betty (54) is presumed to have
#   survived the older Alex (57) under Hong Kong law, making her the sole owner upon his death.
# - Flat B and other assets were solely owned by Betty and are part of her estate.
total_estate_value = (
    flat_a_orchid_apt +
    flat_b_peach_apt +
    bank_accounts +
    nyse_shares +
    personal_items +
    jewellery_collection
)

# --- Step 2: Define and calculate the total value of valid legacies to be paid ---

# Legacy Values from the Will (in HKD)
# - Legacy of Flat C to Cindy: The gift fails due to 'ademption' as the flat was sold. Value is 0.
legacy_for_cindy = 0

# - Legacy to friends: The gift fails for 'uncertainty' as the schedule of friends was never created. Value is 0.
legacy_to_friends = 0

# - Legacy to Wills Lawyers & Co: This is a valid pecuniary legacy.
legacy_to_lawyers = 230_000

# - Legacy to RSPCA: The clause was deleted, so we assume the gift was validly revoked. Value is 0.
legacy_to_rspca = 0

total_legacies_paid = legacy_for_cindy + legacy_to_friends + legacy_to_lawyers + legacy_to_rspca


# --- Step 3: Calculate the residuary estate ---
# Residuary estate = Total Estate - Total Valid Legacies
# (Note: This is before deduction of unspecified debts and administrative expenses)
residuary_estate = total_estate_value - total_legacies_paid


# --- Step 4: Print the final calculation and result ---
print("Calculation of Betty's Residuary Estate")
print("==========================================")
print("\nASSETS (Total Gross Estate):")
print(f"  Flat A (by survivorship) : {flat_a_orchid_apt:>12,}")
print(f"  Flat B (sole ownership)  : {flat_b_peach_apt:>12,}")
print(f"  Bank Accounts            : {bank_accounts:>12,}")
print(f"  NYSE Shares              : {nyse_shares:>12,}")
print(f"  Personal Items           : {personal_items:>12,}")
print(f"  Jewellery                : {jewellery_collection:>12,}")
print(f"                           ------------")
print(f"  Total Estate Value       : {total_estate_value:>12,}")

print("\nLEGACIES (Deductions from Estate):")
print(f"  Gift to Cindy (Failed)   : {legacy_for_cindy:>12,}")
print(f"  Gift to Friends (Failed) : {legacy_to_friends:>12,}")
print(f"  Gift to Lawyers (Valid)  : {legacy_to_lawyers:>12,}")
print(f"  Gift to RSPCA (Revoked)  : {legacy_to_rspca:>12,}")
print(f"                           ------------")
print(f"  Total Legacies Paid      : {total_legacies_paid:>12,}")


print("\nFINAL CALCULATION:")
# Final output showing the equation with each number as requested
print(f"{total_estate_value:,} (Total Estate) - {total_legacies_paid:,} (Total Legacies) = {residuary_estate:,}")
print("==========================================")
print(f"\nThe value of Betty's residuary estate is {residuary_estate:,} HKD.")
