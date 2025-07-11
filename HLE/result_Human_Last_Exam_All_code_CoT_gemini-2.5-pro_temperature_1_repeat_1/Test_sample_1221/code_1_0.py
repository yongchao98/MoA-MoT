# Step 1: Define and calculate the total value of assets in Betty's estate.
# Flat A is included due to the commorientes rule (Betty, being younger, is presumed to have survived Alex).
flat_a = 2000000
flat_b = 4000000
bank_accounts = 50000
shares = 30000
personal_items = 20000
jewellery = 500000

total_assets = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery

print("--- Calculation of Total Estate Assets ---")
print(f"Flat A (by survivorship): {flat_a}")
print(f"Flat B: {flat_b}")
print(f"Bank Accounts: {bank_accounts}")
print(f"Shares: {shares}")
print(f"Personal Items: {personal_items}")
print(f"Jewellery Collection: {jewellery}")
print(f"Total Value of Estate: {flat_a} + {flat_b} + {bank_accounts} + {shares} + {personal_items} + {jewellery} = {total_assets}\n")

# Step 2: Identify the value of valid legacies to be paid from the estate.
# Gift of Flat C: Fails (adeemed).
# Gift to friends (500,000): Fails (void for uncertainty).
# Gift to RSPCA (150,000): Fails (revoked/deleted).
legacy_to_lawyers = 230000  # This is the only valid legacy.

total_valid_legacies = legacy_to_lawyers

print("--- Calculation of Liabilities from Will ---")
print(f"Value of valid legacies to be paid out (to Wills Lawyers & Co): {total_valid_legacies}\n")
print("Note: Other legacies failed due to ademption (Flat C), uncertainty (friends), or revocation (RSPCA).\n")

# Step 3: Calculate the residuary estate.
residuary_estate = total_assets - total_valid_legacies

print("--- Final Residuary Estate Calculation ---")
print(f"The residuary estate is calculated as Total Assets - Valid Legacies.")
print(f"Final Equation: {total_assets} - {total_valid_legacies} = {residuary_estate}")
