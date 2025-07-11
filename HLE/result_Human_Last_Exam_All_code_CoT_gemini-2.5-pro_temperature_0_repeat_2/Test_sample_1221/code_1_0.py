# Step 1: Define and sum the assets in Betty's estate.
# Flat A passes to Betty via survivorship as she is younger than Alex.
flat_a = 2000000
flat_b = 4000000
bank_accounts = 50000
shares = 30000
personal_items = 20000
jewellery = 500000

total_assets = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery

# Step 2: Define the gifts from the will and determine their validity.
# Gift of Flat C fails (adeemed).
gift_flat_c = 0
# Gift to friends fails (uncertainty).
gift_to_friends = 0
# Gift to Wills Lawyers & Co is valid.
gift_to_lawyers = 230000
# Gift to RSPCA fails (revoked by deletion).
gift_to_rspca = 0

total_valid_gifts = gift_flat_c + gift_to_friends + gift_to_lawyers + gift_to_rspca

# Step 3: Calculate the residuary estate.
residuary_estate = total_assets - total_valid_gifts

# Step 4: Print the final calculation.
print("Calculation of Betty's Residuary Estate:")
print("=========================================")
print("Total Assets:")
print(f"  Flat A (Orchid Apt): {flat_a:,}")
print(f"  Flat B (Peach Apt):  {flat_b:,}")
print(f"  Bank Accounts:       {bank_accounts:,}")
print(f"  Shares:              {shares:,}")
print(f"  Personal Items:      {personal_items:,}")
print(f"  Jewellery:           {jewellery:,}")
print("-----------------------------------------")
print(f"Total Estate Value:    {total_assets:,} HKD")
print("\nLess Valid Gifts from Will:")
print(f"  Gift to Wills Lawyers: {gift_to_lawyers:,}")
print("-----------------------------------------")
print(f"Total Valid Gifts:     {total_valid_gifts:,} HKD")
print("\nFinal Calculation:")
print(f"{total_assets:,} (Total Assets) - {total_valid_gifts:,} (Valid Gifts) = {residuary_estate:,} HKD")
print("=========================================")
print(f"The total value of Betty's residuary estate is {residuary_estate:,} HKD.")
