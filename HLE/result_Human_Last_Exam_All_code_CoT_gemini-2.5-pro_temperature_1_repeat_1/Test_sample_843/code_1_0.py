import locale

# Use locale for currency formatting
try:
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
except locale.Error:
    locale.setlocale(locale.LC_ALL, 'C') # Fallback

def format_currency(value):
    return locale.currency(value, grouping=True)

# Accredited Investor Thresholds
INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD = 1000000
INDIVIDUAL_NET_ASSETS_THRESHOLD = 5000000
JOINT_INCOME_THRESHOLD = 300000
ENTITY_NET_ASSETS_THRESHOLD = 5000000

print("Analyzing each option based on Accredited Investor criteria in Ontario...\n")

# --- Option A: Limited Partnership ---
print("--- Analysis of Option A ---")
liam_net_financial_assets = 5400000
jack_assets = 18000000
jack_liabilities = 5000000
ace_net_financial_assets = 25000000

jack_net_assets = jack_assets - jack_liabilities
is_liam_ai = liam_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD
is_jack_ai = jack_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD
is_ace_ai = ace_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD

print(f"Liam's net financial assets: {format_currency(liam_net_financial_assets)}. Is Liam an AI? {is_liam_ai}")
print(f"Jack's net assets: {format_currency(jack_assets)} - {format_currency(jack_liabilities)} = {format_currency(jack_net_assets)}. Is Jack an AI? {is_jack_ai}")
print(f"Ace's net financial assets: {format_currency(ace_net_financial_assets)}. Is Ace an AI? {is_ace_ai}")
if is_liam_ai and is_jack_ai and is_ace_ai:
    print("Conclusion: All limited partners are accredited investors. Therefore, the partnership qualifies as an AI under the look-through provision.\n")

# --- Option B: Individual (Joint Income) ---
print("--- Analysis of Option B ---")
individual_income_2019 = 150000
spouse_income_2019 = 170000
individual_income_2020 = 175000
spouse_income_2020 = 175000

joint_income_2019 = individual_income_2019 + spouse_income_2019
joint_income_2020 = individual_income_2020 + spouse_income_2020
meets_2019_test = joint_income_2019 > JOINT_INCOME_THRESHOLD
meets_2020_test = joint_income_2020 > JOINT_INCOME_THRESHOLD

print(f"2019 Joint Income: {format_currency(individual_income_2019)} + {format_currency(spouse_income_2019)} = {format_currency(joint_income_2019)}. Met threshold? {meets_2019_test}")
print(f"2020 Joint Income: {format_currency(individual_income_2020)} + {format_currency(spouse_income_2020)} = {format_currency(joint_income_2020)}. Met threshold? {meets_2020_test}")
if meets_2019_test and meets_2020_test:
    print("Conclusion: The individual qualifies as an AI under the joint net income test.\n")

# --- Option C: Individual (Net Assets) ---
print("--- Analysis of Option C ---")
joint_assets = 6000000
joint_liabilities = 1000000
joint_net_assets = joint_assets - joint_liabilities
meets_net_asset_test = joint_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD

print(f"Joint Net Assets: {format_currency(joint_assets)} - {format_currency(joint_liabilities)} = {format_currency(joint_net_assets)}. Met threshold? {meets_net_asset_test}")
if meets_net_asset_test:
    print("Conclusion: The individual qualifies as an AI under the net assets test.\n")

# --- Option D: Corporation (Jose and James) ---
print("--- Analysis of Option D ---")
jose_net_financial_assets = 100000000
transfer_pct = 0.10
corp_d_net_assets = jose_net_financial_assets * transfer_pct
meets_corp_asset_test = corp_d_net_assets > ENTITY_NET_ASSETS_THRESHOLD

print(f"Corporation's Net Assets: {format_currency(jose_net_financial_assets)} * {transfer_pct:.0%} = {format_currency(corp_d_net_assets)}. Met threshold? {meets_corp_asset_test}")
if meets_corp_asset_test:
    print("Conclusion: The corporation has net assets exceeding $5,000,000 and qualifies as an AI.\n")

# --- Option E: Corporation (Alex and Bob) ---
print("--- Analysis of Option E ---")
corp_e_assets = 5500000
corp_e_liabilities = 300000
corp_e_net_assets = corp_e_assets - corp_e_liabilities
meets_corp_e_asset_test = corp_e_net_assets > ENTITY_NET_ASSETS_THRESHOLD

print(f"Corporation's Net Assets: {format_currency(corp_e_assets)} - {format_currency(corp_e_liabilities)} = {format_currency(corp_e_net_assets)}. Met threshold? {meets_corp_e_asset_test}")
print("While the corporation meets the net asset test on its face, we must examine its owners.")

alex_nfa = 900000
alex_na = 3000000
is_alex_ai = (alex_nfa > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD) or (alex_na >= INDIVIDUAL_NET_ASSETS_THRESHOLD)
print(f"Alex's Net Financial Assets: {format_currency(alex_nfa)} (Threshold: {format_currency(INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD)}).")
print(f"Alex's Net Assets: {format_currency(alex_na)} (Threshold: {format_currency(INDIVIDUAL_NET_ASSETS_THRESHOLD)}).")
print(f"Is Alex an AI? {is_alex_ai}")

bob_income = 41000
bob_nfa = 75000
# Bob's total assets do not exceed his total liabilities, so net assets <= 0
bob_na = 0
is_bob_ai = (bob_nfa > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD) or (bob_na >= INDIVIDUAL_NET_ASSETS_THRESHOLD)
print(f"Bob's Net Financial Assets: {format_currency(bob_nfa)} (Threshold: {format_currency(INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD)}). Bob's Net Assets <= 0.")
print(f"Is Bob an AI? {is_bob_ai}")

if not is_alex_ai and not is_bob_ai:
    print("\nFinal Conclusion: Both owners, Alex and Bob, are not accredited investors. A corporation formed by non-accredited individuals solely to meet the $5,000,000 net asset test is generally not permitted to be treated as an Accredited Investor. Therefore, this corporation would not qualify.")