# Define the financial details for the individual and spouse from option C.
total_assets = 6000000.00
total_liabilities = 1000000.00

# Define the threshold for the Net Asset Test for an individual (alone or with a spouse).
# The rule requires net assets to EXCEED this value.
net_asset_threshold = 5000000.00

# Calculate the net assets.
net_assets = total_assets - total_liabilities

# Check if the calculated net assets meet the accredited investor criteria.
is_accredited = net_assets > net_asset_threshold

# Print the analysis.
print("--- Accredited Investor Test: Net Asset Rule for an Individual ---")
print(f"Total Assets: ${total_assets:,.2f}")
print(f"Total Liabilities: ${total_liabilities:,.2f}")
print("---------------------------------------------------------------")
# The final code must output each number in the final equation.
print(f"Net Asset Calculation: ${total_assets:,.2f} - ${total_liabilities:,.2f} = ${net_assets:,.2f}")
print("\n--- Compliance Check ---")
print(f"Required Net Assets: Must exceed ${net_asset_threshold:,.2f}")
print(f"Their Net Assets: ${net_assets:,.2f}")
print(f"Do they meet the requirement? {net_assets} > {net_asset_threshold} is {is_accredited}")

if not is_accredited:
    print("\nConclusion: The individual/couple would NOT be classified as an Accredited Investor.")
    print("Reason: Their net assets are exactly $5,000,000.00, which does not exceed the required threshold of $5,000,000.00.")
else:
    print("\nConclusion: The individual/couple would be classified as an Accredited Investor.")
