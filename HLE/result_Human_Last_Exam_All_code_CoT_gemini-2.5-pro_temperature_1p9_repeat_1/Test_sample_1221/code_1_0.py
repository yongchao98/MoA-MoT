# Step 1: Define the values of all assets that form Betty's estate.
# Flat A passed entirely to Betty due to the commorientes rule (younger survives elder).
flat_a = 2000000
# Flat B was owned solely by Betty.
flat_b = 4000000
# Other assets.
cash_in_bank = 50000
shares = 30000
personal_items = 20000
jewellery = 500000

# Create a list of all assets for summation and printing.
assets = [flat_a, flat_b, cash_in_bank, shares, personal_items, jewellery]

# Calculate the total value of the estate.
total_estate = sum(assets)

# Step 2: Define the values of all legacies mentioned in the will.
# Gift of Flat C to Cindy: Fails (ademption) as the flat was sold. Value = 0.
gift_to_cindy = 0
# Gift to friends: Fails (uncertainty) as the schedule was not created. Value = 0.
gift_to_friends = 0
# Gift to RSPCA: Fails as the clause was deleted. Value = 0.
gift_to_rspca = 0
# Gift to Wills Lawyers & Co: This is a valid legacy.
gift_to_lawyers = 230000

# Create a list of only the valid legacies that must be paid out.
valid_legacies_values = [gift_to_lawyers]

# Calculate the total value of valid legacies.
total_valid_legacies = sum(valid_legacies_values)

# Step 3: Calculate the residuary estate.
# Residuary Estate = Total Estate - Total Valid Legacies.
residuary_estate = total_estate - total_valid_legacies

# Print the final calculation showing all the components.
# The equation shows the sum of all assets minus the sum of all valid legacies.
assets_str = ' + '.join(map(str, assets))
legacies_str = ' + '.join(map(str, valid_legacies_values))

print(f"Calculation of Betty's Residuary Estate (in HKD):")
print(f"({assets_str}) - ({legacies_str}) = {residuary_estate}")
