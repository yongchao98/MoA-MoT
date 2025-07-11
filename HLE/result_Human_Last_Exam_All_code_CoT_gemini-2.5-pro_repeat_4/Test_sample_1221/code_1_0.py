# Step 1: Define the value of all assets in Betty's estate.

# Flat A is included in full due to the commorientes rule (the younger, Betty, is presumed to survive the elder, Alex)
# which makes her the sole owner by survivorship at the time of her death.
flat_a = 2000000
# Flat B was owned solely by Betty.
flat_b = 4000000
# Other personal assets.
bank_accounts = 50000
shares = 30000
personal_items = 20000
jewellery = 500000

# Calculate the total value of the gross estate.
total_assets = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery

# Step 2: Define the value of the valid legacies to be paid out from the estate.

# The gift of Flat C to Cindy is adeemed (it fails) because Betty did not own it at death.
gift_to_cindy = 0
# The gift to friends fails for uncertainty as the schedule was never created.
legacy_to_friends = 0
# The gift to Wills Lawyers & Co is a valid pecuniary legacy.
legacy_to_lawyers = 230000
# The gift to RSPCA was revoked by deletion.
legacy_to_rspca = 0

# Calculate the total value of legacies that must be paid.
total_payouts = gift_to_cindy + legacy_to_friends + legacy_to_lawyers + legacy_to_rspca

# Step 3: Calculate the residuary estate.
# The residuary estate is what remains after all valid legacies and debts are paid.
# We assume no other debts for this calculation as none were specified.
residuary_estate = total_assets - total_payouts

# Step 4: Print the final calculation, showing all the numbers.
print("Calculation of Betty's Residuary Estate:")
print("=" * 40)

print("\n--- Total Assets ---")
print(f"Value of Flat A (by survivorship): {flat_a:,}")
print(f"Value of Flat B: {flat_b:,}")
print(f"Value of Bank Accounts: {bank_accounts:,}")
print(f"Value of Shares: {shares:,}")
print(f"Value of Personal Items: {personal_items:,}")
print(f"Value of Jewellery: {jewellery:,}")
print("--------------------")
print(f"Equation for Total Assets: {flat_a} + {flat_b} + {bank_accounts} + {shares} + {personal_items} + {jewellery}")
print(f"Total Assets = {total_assets:,} HKD")

print("\n--- Valid Payouts (Legacies) ---")
print(f"Legacy to Wills Lawyers & Co: {legacy_to_lawyers:,}")
print("--------------------")
print(f"Total Valid Payouts = {total_payouts:,} HKD")
print("(Note: Other legacies failed due to ademption, uncertainty, or revocation.)")

print("\n--- Final Residuary Estate Calculation ---")
print("Equation: Total Assets - Total Valid Payouts")
print(f"Final Equation: {total_assets} - {total_payouts} = {residuary_estate}")
print("--------------------")
print(f"Betty's Residuary Estate = {residuary_estate:,} HKD")
