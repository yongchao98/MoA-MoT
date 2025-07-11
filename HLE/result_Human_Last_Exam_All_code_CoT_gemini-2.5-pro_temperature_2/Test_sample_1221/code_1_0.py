# --- Step 1: Calculate the total value of Betty's assets ---

# Asset values in HKD
flat_a_orchid = 2000000 # Passed entirely to Betty via survivorship
flat_b_peach = 4000000
bank_accounts = 50000
shares_nyse = 30000
personal_items = 20000
jewellery = 500000

# Summing up all assets
total_assets = flat_a_orchid + flat_b_peach + bank_accounts + shares_nyse + personal_items + jewellery

# --- Step 2: Identify valid legacies to be paid out ---

# Clause 4 (Flat C): Gift fails (ademption)
gift_cindy = 0

# Clause 5a (Friends): Gift fails (uncertainty)
gift_friends = 0

# Clause 5b (Wills Lawyers & Co): Gift is valid
gift_lawyers = 230000

# Clause 6 (RSPCA): Gift fails (revoked)
gift_rspca = 0

# Sum of all valid payouts (legacies) from the estate
total_payouts = gift_cindy + gift_friends + gift_lawyers + gift_rspca

# --- Step 3: Calculate the residuary estate ---

# Residuary estate is total assets minus total payouts
# Note: Debts and administrative/funeral costs are not specified and thus omitted from this calculation.
residuary_estate = total_assets - total_payouts

# --- Step 4: Print the final calculation and result ---

print("Calculation of Betty's Residuary Estate:")
print("Total Assets:")
print(f"  Flat A (survivorship): {flat_a_orchid}")
print(f"  Flat B: {flat_b_peach}")
print(f"  Bank Accounts: {bank_accounts}")
print(f"  Shares: {shares_nyse}")
print(f"  Personal Items: {personal_items}")
print(f"  Jewellery: {jewellery}")
print(f"-------------------------------")
print(f"Total Assets = {total_assets}\n")

print("Valid Legacies to be Paid Out:")
print(f"  Legacy to Wills Lawyers & Co: {gift_lawyers}")
print(f"-------------------------------")
print(f"Total Payouts = {total_payouts}\n")

print("Final Residuary Estate Calculation:")
print(f"{total_assets} (Total Assets) - {total_payouts} (Total Payouts) = {residuary_estate}")

print(f"\nBetty's residuary estate is {residuary_estate} HKD.")
<<<6370000>>>