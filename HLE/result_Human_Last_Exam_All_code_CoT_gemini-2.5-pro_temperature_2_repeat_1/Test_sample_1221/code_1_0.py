# Plan:
# 1. Define all of Betty's assets at the time of her death.
#    - Flat A is included fully due to the rule of survivorship in joint tenancy, as Betty is presumed to have survived Alex.
#    - Flat B, cash, shares, personal items, and jewellery are all part of her estate.
# 2. Define all valid pecuniary legacies (gifts of money) from her will that must be paid out.
#    - The gift of Flat C fails (adeemed).
#    - The gift to friends fails (uncertainty of object).
#    - The gift to RSPCA fails (revoked/deleted).
#    - The gift to Wills Lawyers & Co is valid.
# 3. Calculate total assets.
# 4. Calculate total liabilities (valid gifts).
# 5. Calculate the residuary estate by subtracting total liabilities from total assets.
# 6. Print the full equation showing all numbers.

# Assets (in HKD)
# Flat A (Orchid Apartment): Value passes to Betty's estate via survivorship.
flat_a = 2000000
# Flat B (Peach Apartment): Owned solely by Betty.
flat_b = 4000000
# Cash in bank accounts.
cash_in_bank = 50000
# NYSE Shares.
shares = 30000
# Personal items.
personal_items = 20000
# Jewellery collection.
jewellery = 500000

# Calculate total assets
total_assets = flat_a + flat_b + cash_in_bank + shares + personal_items + jewellery

# Valid Liabilities from the Will (in HKD)
# The only valid legacy to be paid out.
gift_to_lawyers = 230000

# Total liabilities from the will.
total_liabilities = gift_to_lawyers

# Calculate the residuary estate
residuary_estate = total_assets - total_liabilities

# Print the final equation as requested
print("Calculation of Betty's residuary estate:")
print("Assets - Liabilities = Residuary Estate")
print(f"({flat_a} + {flat_b} + {cash_in_bank} + {shares} + {personal_items} + {jewellery}) - {gift_to_lawyers} = {residuary_estate}")
print(f"Total: {residuary_estate} HKD")
<<<6370000>>>