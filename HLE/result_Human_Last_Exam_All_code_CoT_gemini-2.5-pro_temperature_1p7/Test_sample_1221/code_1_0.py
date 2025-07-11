# Step 1: Define the value of assets in Betty's estate
# Flat A: Passed to Betty via survivorship as she was younger than Alex.
flat_a = 2000000
# Flat B: Owned solely by Betty.
flat_b = 4000000
# Cash in bank accounts.
bank_accounts = 50000
# Value of shares.
shares = 30000
# Value of personal items.
personal_items = 20000
# Value of jewellery collection.
jewellery = 500000

# Calculate the total value of the gross estate
total_gross_estate = flat_a + flat_b + bank_accounts + shares + personal_items + jewellery
print("Calculating the Total Gross Estate:")
print(f"Flat A (via survivorship): {flat_a:,} HKD")
print(f"Flat B (sole ownership): {flat_b:,} HKD")
print(f"Bank Accounts: {bank_accounts:,} HKD")
print(f"Shares: {shares:,} HKD")
print(f"Personal Items: {personal_items:,} HKD")
print(f"Jewellery: {jewellery:,} HKD")
print(f"Total Gross Estate = {total_gross_estate:,} HKD\n")


# Step 2: Define the value of legacies to be paid out from the estate
# Gift of Flat C failed due to ademption (asset was sold).
gift_to_cindy = 0
# Gift to friends failed due to uncertainty of objects (beneficiaries not named).
gift_to_friends = 0
# Gift to Wills Lawyers & Co is a valid pecuniary legacy.
gift_to_lawyers = 230000
# Gift to RSPCA was revoked (clause deleted).
gift_to_rspca = 0

# Calculate the total value of legacies to be deducted
total_legacies = gift_to_cindy + gift_to_friends + gift_to_lawyers + gift_to_rspca
print("Calculating Total Deductions from Valid Gifts:")
print(f"Gift of Flat C (failed): {gift_to_cindy:,} HKD")
print(f"Gift to Friends (failed): {gift_to_friends:,} HKD")
print(f"Gift to Wills Lawyers & Co (valid): {gift_to_lawyers:,} HKD")
print(f"Gift to RSPCA (revoked): {gift_to_rspca:,} HKD")
print(f"Total Deductions = {total_legacies:,} HKD\n")


# Step 3: Calculate the residuary estate
# Residuary estate = Total Gross Estate - Total Legacies.
# Note: This is the value before deducting unspecified debts and administrative costs.
residuary_estate = total_gross_estate - total_legacies

print("Final Calculation of the Residuary Estate:")
print(f"Residuary Estate = Total Gross Estate - Total Deductions")
print(f"Residuary Estate = {total_gross_estate:,} - {total_legacies:,}")
print(f"The total value of Betty's residuary estate is: {residuary_estate:,} HKD")
<<<6370000>>>