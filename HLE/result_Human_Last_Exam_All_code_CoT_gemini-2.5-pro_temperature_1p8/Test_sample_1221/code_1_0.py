# --- Asset Calculation ---

# Flat A (Orchid Apartment): 2,000,000 HKD
# Under the rule of survivorship and the presumption of death (younger survives elder),
# the full value of the jointly owned flat falls into Betty's estate.
flat_a = 2000000

# Flat B (Peach Apartment): 4,000,000 HKD
# This flat was solely owned by Betty.
flat_b = 4000000

# Bank Accounts: 50,000 HKD
bank_accounts = 50000

# NYSE Shares: 30,000 HKD
shares = 30000

# Personal Items: 20,000 HKD
personal_items = 20000

# Jewellery Collection: 500,000 HKD
jewellery = 500000

# --- Legacy (Payout) Calculation ---

# Gift of Flat C to Cindy: 0 HKD
# The gift is 'adeemed' (fails) because Betty no longer owned the property at her death.
gift_cindy = 0

# Gift to Friends: 0 HKD
# This gift fails for 'uncertainty of objects' as the schedule of friends was never created.
gift_friends = 0

# Gift to Wills Lawyers & Co: 230,000 HKD
# This is a valid pecuniary legacy.
gift_lawyers = 230000

# Gift to RSPCA: 150,000 HKD
# A simple deletion of a clause is not legally effective without proper formalities.
# Therefore, the legacy is presumed to be valid.
gift_rspca = 150000

# --- Final Calculation ---

# Total value of all assets in the estate
total_assets = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery

# Total value of all valid legacies to be paid out
total_payouts = gift_cindy + gift_friends + gift_lawyers + gift_rspca

# The residuary estate is the total assets minus the total payouts.
residuary_estate = total_assets - total_payouts

# Print the equation and the final result
print(f"Total Assets ({flat_a:,} + {flat_b:,} + {bank_accounts:,} + {shares:,} + {personal_items:,} + {jewellery:,}) = {total_assets:,} HKD")
print(f"Total Valid Legacies ({gift_cindy:,} + {gift_friends:,} + {gift_lawyers:,} + {gift_rspca:,}) = {total_payouts:,} HKD")
print(f"Residuary Estate Calculation: {total_assets:,} - {total_payouts:,} = {residuary_estate:,} HKD")
print(f"\nThe total value of Betty's residuary estate is {residuary_estate:,} HKD.")
<<<6220000>>>