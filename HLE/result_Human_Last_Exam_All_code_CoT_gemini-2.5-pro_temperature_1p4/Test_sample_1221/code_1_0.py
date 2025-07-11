import locale

# --- Betty's Assets ---

# Flat A is included in the estate due to the commorientes rule (Betty, being younger, is presumed to have survived Alex).
flat_a_value = 2000000
# Flat B was solely owned by Betty.
flat_b_value = 4000000
# Other personal assets.
cash_value = 50000
shares_value = 30000
personal_items_value = 20000
jewellery_value = 500000

# --- Will Analysis for Payouts ---

# Legacy to Cindy (Flat C): Failed (adeemed) because the flat was sold. Payout = 0.
payout_cindy = 0
# Legacy to Friends: Failed due to uncertainty (schedule not created). Payout = 0.
payout_friends = 0
# Legacy to RSPCA: Failed because the clause was deleted. Payout = 0.
payout_rspca = 0
# Legacy to Wills Lawyers & Co: This is a valid legacy.
payout_lawyers = 230000

# --- Calculation ---

# 1. Sum up all assets that form the estate.
total_gross_estate = flat_a_value + flat_b_value + cash_value + shares_value + personal_items_value + jewellery_value

# 2. Sum up all valid payouts from the estate.
total_payouts = payout_cindy + payout_friends + payout_rspca + payout_lawyers

# 3. Calculate the residuary estate.
residuary_estate = total_gross_estate - total_payouts

# --- Print the results ---

# Use locale for better number formatting (e.g., adding commas)
try:
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
except locale.Error:
    try:
        locale.setlocale(locale.LC_ALL, 'English_United States.1252')
    except locale.Error:
        pass # Fallback to default formatting if locales are not available

print("Calculating the Residuary Estate of Betty:")
print("-" * 50)
print("Total Gross Estate Calculation:")
print(f"  Value of Flat A (Orchid Apmt): {locale.format_string('%d', flat_a_value, grouping=True)} HKD")
print(f"  Value of Flat B (Peach Apmt):  {locale.format_string('%d', flat_b_value, grouping=True)} HKD")
print(f"  Value of Cash and Shares:     {locale.format_string('%d', (cash_value + shares_value), grouping=True)} HKD")
print(f"  Value of Items and Jewellery: {locale.format_string('%d', (personal_items_value + jewellery_value), grouping=True)} HKD")
print(f"Total Gross Estate:               {locale.format_string('%d', total_gross_estate, grouping=True)} HKD")
print("-" * 50)
print("Total Valid Payouts (Legacies):")
print(f"  Valid Legacy to Wills Lawyers & Co: {locale.format_string('%d', payout_lawyers, grouping=True)} HKD")
print(f"Total Payouts:                      {locale.format_string('%d', total_payouts, grouping=True)} HKD")
print("-" * 50)
print("Final Residuary Estate Calculation:")
print(f"  {locale.format_string('%d', total_gross_estate, grouping=True)} (Total Estate) - {locale.format_string('%d', total_payouts, grouping=True)} (Total Payouts)")
print(f"Result: {locale.format_string('%d', residuary_estate, grouping=True)} HKD\n")
print("Final Equation:")
print(f"{flat_a_value} + {flat_b_value} + {cash_value} + {shares_value} + {personal_items_value} + {jewellery_value} - {payout_lawyers} = {residuary_estate}")

<<<6370000>>>