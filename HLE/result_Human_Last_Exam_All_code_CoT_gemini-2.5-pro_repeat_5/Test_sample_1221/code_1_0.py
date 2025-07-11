# Step 1: Define the value of Betty's assets.
# Flat A is included due to the commorientes rule, where Betty (younger) is presumed to have survived Alex (elder).
flat_a_value = 2000000
flat_b_value = 4000000
cash_in_bank = 50000
shares_value = 30000
personal_items_value = 20000
jewellery_value = 500000

# Calculate the total value of the estate.
total_estate_assets = [
    flat_a_value,
    flat_b_value,
    cash_in_bank,
    shares_value,
    personal_items_value,
    jewellery_value
]
total_estate_value = sum(total_estate_assets)

# Step 2: Define the value of gifts from the will and assess their validity.
# Gift of Flat C to Cindy: Fails (Ademption)
gift_to_cindy = 0
# Gift to friends in schedule: Fails (Uncertainty of Objects)
gift_to_friends = 0
# Gift to Wills Lawyers & Co: Valid
gift_to_lawyers = 230000
# Gift to RSPCA: Fails (Revocation by deletion)
gift_to_rspca = 0

# Calculate the total value of valid gifts to be paid out.
total_valid_gifts = gift_to_cindy + gift_to_friends + gift_to_lawyers + gift_to_rspca

# Step 3: Calculate the residuary estate.
# Residuary Estate = Total Estate Value - Total Valid Gifts
residuary_estate = total_estate_value - total_valid_gifts

# Output the results and the final equation.
print("Calculating Betty's Residuary Estate:")
print("-" * 40)
print("Assets:")
print(f"  - Flat A (Orchid Apt): {flat_a_value:,} HKD")
print(f"  - Flat B (Peach Apt):  {flat_b_value:,} HKD")
print(f"  - Bank Accounts:       {cash_in_bank:,} HKD")
print(f"  - NYSE Shares:         {shares_value:,} HKD")
print(f"  - Personal Items:      {personal_items_value:,} HKD")
print(f"  - Jewellery:           {jewellery_value:,} HKD")
print(f"Total Estate Value: {total_estate_value:,} HKD\n")

print("Gifts from Will:")
print(f"  - Gift to Cindy (Flat C): FAILED")
print(f"  - Gift to Friends (Schedule): FAILED")
print(f"  - Gift to Wills Lawyers: {gift_to_lawyers:,} HKD")
print(f"  - Gift to RSPCA: FAILED")
print(f"Total Value of Valid Gifts: {total_valid_gifts:,} HKD\n")

print("-" * 40)
print("Final Calculation of Residuary Estate:")
print(f"Total Estate ({total_estate_value:,}) - Total Valid Gifts ({total_valid_gifts:,})")
print(f"= {residuary_estate:,} HKD")
print("-" * 40)
print(f"Betty's residuary estate is {residuary_estate:,} HKD.")