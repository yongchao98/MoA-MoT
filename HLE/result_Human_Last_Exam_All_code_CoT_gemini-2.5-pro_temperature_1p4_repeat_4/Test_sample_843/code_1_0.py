import sys

# Suppress scientific notation for cleaner output
def format_currency(value):
    return f"${value:,.2f}"

# Accredited Investor (AI) Thresholds
financial_asset_threshold = 1000000
net_asset_threshold = 5000000
individual_income_threshold = 200000
joint_income_threshold = 300000
corporation_net_asset_threshold = 5000000

print("Analysis of each option based on Accredited Investor rules in Ontario:")
print("-" * 60)

# --- Option A: The Limited Partnership ---
print("Option A: The Limited Partnership")
liam_net_financial_assets = 5400000.00
jack_total_assets = 18000000.00
jack_total_liabilities = 5000000.00
ace_net_financial_assets = 25000000.00
gift_per_partner_to_gp = 2000000.00 # Note: The prompt has a typo "2,000,0000.00", assumed to be 2,000,000

# Check Liam
liam_is_ai = liam_net_financial_assets > financial_asset_threshold
print(f"Liam's net financial assets of {format_currency(liam_net_financial_assets)} are greater than the {format_currency(financial_asset_threshold)} threshold. Liam is an AI: {liam_is_ai}")

# Check Jack
jack_net_assets = jack_total_assets - jack_total_liabilities
jack_is_ai = jack_net_assets >= net_asset_threshold
print(f"Jack's net assets ({format_currency(jack_total_assets)} - {format_currency(jack_total_liabilities)}) are {format_currency(jack_net_assets)}, which meets the {format_currency(net_asset_threshold)} threshold. Jack is an AI: {jack_is_ai}")

# Check Ace
ace_is_ai = ace_net_financial_assets > financial_asset_threshold
print(f"Ace's net financial assets of {format_currency(ace_net_financial_assets)} are greater than the {format_currency(financial_asset_threshold)} threshold. Ace is an AI: {ace_is_ai}")

# Check General Partner (GP)
gp_net_assets = 3 * gift_per_partner_to_gp
gp_is_ai = gp_net_assets >= corporation_net_asset_threshold
print(f"The GP's net assets are 3 * {format_currency(gift_per_partner_to_gp)} = {format_currency(gp_net_assets)}, which meets the {format_currency(corporation_net_asset_threshold)} threshold. The GP is an AI: {gp_is_ai}")

lp_is_ai = liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai
print(f"Conclusion: Since all owners are accredited investors, the Limited Partnership IS an Accredited Investor based on the look-through provision.")
print("-" * 60)

# --- Option B: The Individual (Joint Income) ---
print("Option B: The Individual (Joint Income)")
individual_income_2019 = 150000.00
spouse_income_2019 = 170000.00
individual_income_2020 = 175000.00
spouse_income_2020 = 175000.00

joint_income_2019 = individual_income_2019 + spouse_income_2019
joint_income_2020 = individual_income_2020 + spouse_income_2020
passes_2019 = joint_income_2019 > joint_income_threshold
passes_2020 = joint_income_2020 > joint_income_threshold

print(f"2019 Joint Income: {format_currency(individual_income_2019)} + {format_currency(spouse_income_2019)} = {format_currency(joint_income_2019)}. Passes (> {format_currency(joint_income_threshold)}): {passes_2019}")
print(f"2020 Joint Income: {format_currency(individual_income_2020)} + {format_currency(spouse_income_2020)} = {format_currency(joint_income_2020)}. Passes (> {format_currency(joint_income_threshold)}): {passes_2020}")
individual_b_is_ai = passes_2019 and passes_2020
print(f"Conclusion: The individual meets the joint income test for the two most recent years. The individual IS an Accredited Investor.")
print("-" * 60)

# --- Option C: The Individual (Joint Net Assets) ---
print("Option C: The Individual (Joint Net Assets)")
total_assets = 6000000.00
total_liabilities = 1000000.00
net_assets = total_assets - total_liabilities
individual_c_is_ai = net_assets >= net_asset_threshold
print(f"Joint Net Assets: {format_currency(total_assets)} - {format_currency(total_liabilities)} = {format_currency(net_assets)}.")
print(f"This meets the 'at least {format_currency(net_asset_threshold)}' test. The individual IS an Accredited Investor.")
print("-" * 60)

# --- Option D: The Corporation (Jose and James) ---
print("Option D: The Corporation (Jose and James)")
# Check Owners
jose_net_financial_assets = 100000000.00
james_joint_income = 75000.00 + 210000.00
print(f"Jose is an AI (net financial assets of {format_currency(jose_net_financial_assets)}). James is not an AI (joint income of {format_currency(james_joint_income)} is less than {format_currency(joint_income_threshold)}).")
print("The corporation fails the look-through test.")
# Check Corporation's Net Assets
transfer_percent = 0.10
corp_d_net_assets = jose_net_financial_assets * transfer_percent
corp_d_is_ai = corp_d_net_assets >= corporation_net_asset_threshold
print(f"The corporation's net assets are 10% of {format_currency(jose_net_financial_assets)} = {format_currency(corp_d_net_assets)}.")
print(f"This meets the {format_currency(corporation_net_asset_threshold)} threshold. The corporation IS an Accredited Investor on its own.")
print("-" * 60)

# --- Option E: The Corporation (Alex and Bob) ---
print("Option E: The Corporation (Alex and Bob)")
# Check Owners
alex_nfa = 900000.00
alex_na = 3000000.00
alex_is_ai = (alex_nfa > financial_asset_threshold) or (alex_na >= net_asset_threshold)
print(f"Alex is not an AI (fails financial asset and net asset tests).")
bob_is_ai = False # Fails all income and asset tests based on description.
print(f"Bob is not an AI (fails income and asset tests).")
print("The corporation fails the look-through test since both owners are not AIs.")

# Check Corporation's Net Assets
corp_e_gic = 5500000.00
corp_e_liabilities = 300000.00
corp_e_min_net_assets = corp_e_gic - corp_e_liabilities
print(f"The corporation's net assets are at least {format_currency(corp_e_gic)} - {format_currency(corp_e_liabilities)} = {format_currency(corp_e_min_net_assets)}.")
print(f"Technically, this value is greater than the {format_currency(corporation_net_asset_threshold)} threshold.")
print("\nFinal Conclusion:")
print("Although the corporation in Option E technically meets the net asset test, it is a vehicle created by two non-accredited investors to pool assets.")
print("Under securities law anti-avoidance principles (Companion Policy 45-106CP), such a structure, created solely to meet the threshold, would NOT be classified as an Accredited Investor.")