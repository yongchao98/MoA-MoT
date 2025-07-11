import locale

# Use locale to format numbers with commas for readability
locale.setlocale(locale.LC_ALL, '')

# Define accredited investor thresholds
NFA_THRESHOLD_INDIVIDUAL = 1_000_000
NET_ASSET_THRESHOLD_INDIVIDUAL = 5_000_000
INCOME_THRESHOLD_INDIVIDUAL = 200_000
INCOME_THRESHOLD_SPOUSAL = 300_000
NET_ASSET_THRESHOLD_ENTITY = 5_000_000

print("Analyzing each option based on Accredited Investor definitions in Ontario.\n")

# --- Analysis of Option A ---
print("--- Option A: Limited Partnership ---")
liam_nfa = 5_400_000
jack_net_assets = 18_000_000 - 5_000_000
ace_nfa = 25_000_000
# The general partner is a corporation owned by Liam, Jack, and Ace.
# Since all three individuals qualify as AIs on their own, their corporation also qualifies.
# Because all limited partners and the general partner are accredited investors, the partnership itself qualifies.
print(f"Liam's net financial assets: ${locale.format_string('%d', liam_nfa, grouping=True)} > ${locale.format_string('%d', NFA_THRESHOLD_INDIVIDUAL, grouping=True)} (Qualifies)")
print(f"Jack's net assets: ${locale.format_string('%d', 18_000_000, grouping=True)} - ${locale.format_string('%d', 5_000_000, grouping=True)} = ${locale.format_string('%d', jack_net_assets, grouping=True)} >= ${locale.format_string('%d', NET_ASSET_THRESHOLD_INDIVIDUAL, grouping=True)} (Qualifies)")
print(f"Ace's net financial assets: ${locale.format_string('%d', ace_nfa, grouping=True)} > ${locale.format_string('%d', NFA_THRESHOLD_INDIVIDUAL, grouping=True)} (Qualifies)")
print("Result for A: All owners are AIs, so the partnership IS an Accredited Investor.\n")

# --- Analysis of Option B ---
print("--- Option B: Individual (Spousal Income Test) ---")
combined_income_2019 = 150_000 + 170_000
combined_income_2020 = 175_000 + 175_000
print(f"2019 Combined Income: ${locale.format_string('%d', 150000, grouping=True)} + ${locale.format_string('%d', 170000, grouping=True)} = ${locale.format_string('%d', combined_income_2019, grouping=True)} > ${locale.format_string('%d', INCOME_THRESHOLD_SPOUSAL, grouping=True)}")
print(f"2020 Combined Income: ${locale.format_string('%d', 175000, grouping=True)} + ${locale.format_string('%d', 175000, grouping=True)} = ${locale.format_string('%d', combined_income_2020, grouping=True)} > ${locale.format_string('%d', INCOME_THRESHOLD_SPOUSAL, grouping=True)}")
print("Result for B: Meets the spousal income test. IS an Accredited Investor.\n")

# --- Analysis of Option C ---
print("--- Option C: Individual (Net Asset Test) ---")
individual_c_net_assets = 6_000_000 - 1_000_000
print(f"Net Assets: ${locale.format_string('%d', 6000000, grouping=True)} - ${locale.format_string('%d', 1000000, grouping=True)} = ${locale.format_string('%d', individual_c_net_assets, grouping=True)}")
print(f"The result is equal to the ${locale.format_string('%d', NET_ASSET_THRESHOLD_INDIVIDUAL, grouping=True)} threshold.")
print("Result for C: Meets the 'at least $5,000,000' net asset test. IS an Accredited Investor.\n")

# --- Analysis of Option D ---
print("--- Option D: Corporation ---")
# Check the owners: Jose is an AI ($100M NFA). James is not.
james_combined_income = 75000 + 210000
print(f"James' combined spousal income: ${locale.format_string('%d', 75000, grouping=True)} + ${locale.format_string('%d', 210000, grouping=True)} = ${locale.format_string('%d', james_combined_income, grouping=True)}, which is less than the ${locale.format_string('%d', INCOME_THRESHOLD_SPOUSAL, grouping=True)} threshold.")
print("Since not all owners are AIs, the corporation fails that test.")
print("However, checking the corporation's own net assets:")
corp_d_net_assets = 100_000_000 * 0.10
print(f"Corporation's Net Assets: 10% of ${locale.format_string('%d', 100000000, grouping=True)} = ${locale.format_string('%d', corp_d_net_assets, grouping=True)}")
print(f"This is greater than the ${locale.format_string('%d', NET_ASSET_THRESHOLD_ENTITY, grouping=True)} entity threshold.")
print("Result for D: Meets the entity net asset test. IS an Accredited Investor.\n")

# --- Analysis of Option E ---
print("--- Option E: Corporation ---")
# Check the owners:
print("Owner Alex: Fails AI tests (Net financial assets of $900k < $1M; Net assets of $3M < $5M).")
print("Owner Bob: Fails AI tests (Low income, negative or zero net worth).")
print("Since both owners are not AIs, the corporation fails the 'all owners are AI' test.")
print("Checking the corporation's own net assets:")
corp_e_net_assets = 5_500_000 - 300_000
print(f"Corporation's Net Assets: ${locale.format_string('%d', 5500000, grouping=True)} - ${locale.format_string('%d', 300000, grouping=True)} = ${locale.format_string('%d', corp_e_net_assets, grouping=True)}")
print(f"Technically, ${locale.format_string('%d', corp_e_net_assets, grouping=True)} is greater than the ${locale.format_string('%d', NET_ASSET_THRESHOLD_ENTITY, grouping=True)} threshold.")
print("However, securities regulations contain anti-avoidance provisions. Since this corporation's only purpose appears to be to pool funds from non-accredited investors to technically meet the AI threshold, a regulator would likely deem this an artificial creation to circumvent the rules.")
print("Result for E: Due to the anti-avoidance provision, this corporation would likely NOT be classified as an Accredited Investor.\n")

final_answer = "E"
print("Final Conclusion: The entity that would not be classified as an Accredited Investor is the one in Option E.")
print(f'<<<{final_answer}>>>')
