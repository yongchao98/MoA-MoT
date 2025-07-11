import locale

# Use locale to format numbers with commas for readability
locale.setlocale(locale.LC_ALL, '')

# Define Accredited Investor thresholds
IND_NFA_THRESHOLD = 1_000_000
IND_NA_THRESHOLD = 5_000_000
IND_INCOME_THRESHOLD = 200_000
JOINT_INCOME_THRESHOLD = 300_000
COMPANY_NA_THRESHOLD = 5_000_000

print("Analysis of Accredited Investor Status as of Jan 2021 (Ontario)")
print("-" * 70)

# --- Option A: Limited Partnership ---
print("A. Limited Partnership:")
liam_nfa = 5_400_000
jack_assets = 18_000_000
jack_liabilities = 5_000_000
jack_na = jack_assets - jack_liabilities
ace_nfa = 25_000_000
# The GP is owned by Liam, Jack, and Ace, who are all accredited investors, so the GP is also accredited.
all_owners_are_ai = (liam_nfa > IND_NFA_THRESHOLD) and (jack_na >= IND_NA_THRESHOLD) and (ace_nfa > IND_NFA_THRESHOLD)
print(f"  - Liam's Net Financial Assets: ${locale.format_string('%d', liam_nfa, True)} > ${locale.format_string('%d', IND_NFA_THRESHOLD, True)} -> Qualifies")
print(f"  - Jack's Net Assets: ${locale.format_string('%d', jack_assets, True)} - ${locale.format_string('%d', jack_liabilities, True)} = ${locale.format_string('%d', jack_na, True)} >= ${locale.format_string('%d', IND_NA_THRESHOLD, True)} -> Qualifies")
print(f"  - Ace's Net Financial Assets: ${locale.format_string('%d', ace_nfa, True)} > ${locale.format_string('%d', IND_NFA_THRESHOLD, True)} -> Qualifies")
print(f"  - General Partner is owned by accredited investors -> Qualifies")
print("--> Result: Option A qualifies as an Accredited Investor.\n")

# --- Option B: Individual (Joint Income) ---
print("B. Individual (Joint Income Test):")
ind_income_2019 = 150_000
spouse_income_2019 = 170_000
joint_income_2019 = ind_income_2019 + spouse_income_2019
ind_income_2020 = 175_000
spouse_income_2020 = 175_000
joint_income_2020 = ind_income_2020 + spouse_income_2020
pass_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
pass_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD
print(f"  - Joint Income 2019: ${locale.format_string('%d', ind_income_2019, True)} + ${locale.format_string('%d', spouse_income_2019, True)} = ${locale.format_string('%d', joint_income_2019, True)} > ${locale.format_string('%d', JOINT_INCOME_THRESHOLD, True)} -> {pass_2019}")
print(f"  - Joint Income 2020: ${locale.format_string('%d', ind_income_2020, True)} + ${locale.format_string('%d', spouse_income_2020, True)} = ${locale.format_string('%d', joint_income_2020, True)} > ${locale.format_string('%d', JOINT_INCOME_THRESHOLD, True)} -> {pass_2020}")
print("--> Result: Option B qualifies as an Accredited Investor.\n")

# --- Option C: Individual (Joint Net Assets) ---
print("C. Individual (Joint Net Asset Test):")
joint_assets = 6_000_000
joint_liabilities = 1_000_000
joint_net_assets = joint_assets - joint_liabilities
passes_na_test = joint_net_assets >= IND_NA_THRESHOLD
print(f"  - Equation: Total Assets ${locale.format_string('%d', joint_assets, True)} - Total Liabilities ${locale.format_string('%d', joint_liabilities, True)} = Net Assets ${locale.format_string('%d', joint_net_assets, True)}")
print(f"  - Test: Is ${locale.format_string('%d', joint_net_assets, True)} at least ${locale.format_string('%d', IND_NA_THRESHOLD, True)}? -> {passes_na_test}")
print("--> Result: Option C qualifies as an Accredited Investor because its net assets are exactly equal to the threshold, meeting the 'at least' condition. However, this is the only option that doesn't clearly exceed its threshold.\n")

# --- Option D: Corporation ---
print("D. Corporation (Jose and James):")
jose_nfa = 100_000_000
corp_d_assets = 0.10 * jose_nfa
passes_corp_na_test = corp_d_assets >= COMPANY_NA_THRESHOLD
print(f"  - Equation: 10% of ${locale.format_string('%d', jose_nfa, True)} = Net Assets ${locale.format_string('%d', corp_d_assets, True)}")
print(f"  - Test: Is ${locale.format_string('%d', corp_d_assets, True)} at least ${locale.format_string('%d', COMPANY_NA_THRESHOLD, True)}? -> {passes_corp_na_test}")
print("--> Result: Option D qualifies as an Accredited Investor.\n")

# --- Option E: Corporation ---
print("E. Corporation (Alex and Bob):")
corp_e_assets = 5_500_000
corp_e_liabilities = 300_000
corp_e_net_assets = corp_e_assets - corp_e_liabilities
passes_corp_e_na_test = corp_e_net_assets >= COMPANY_NA_THRESHOLD
print(f"  - Equation: GIC ${locale.format_string('%d', corp_e_assets, True)} - Liabilities ${locale.format_string('%d', corp_e_liabilities, True)} = Net Assets ${locale.format_string('%d', corp_e_net_assets, True)}")
print(f"  - Test: Is ${locale.format_string('%d', corp_e_net_assets, True)} at least ${locale.format_string('%d', COMPANY_NA_THRESHOLD, True)}? -> {passes_corp_e_na_test}")
print("--> Result: Option E qualifies as an Accredited Investor.\n")
