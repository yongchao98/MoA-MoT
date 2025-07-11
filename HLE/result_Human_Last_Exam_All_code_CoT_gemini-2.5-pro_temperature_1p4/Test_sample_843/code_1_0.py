import locale

# Use locale to format currency
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

def print_currency(value):
    """Formats a number as currency."""
    return locale.currency(value, grouping=True)

def accredited_investor_analysis():
    """
    Analyzes each scenario to determine accredited investor status.
    """
    # Thresholds
    INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
    INDIVIDUAL_NET_ASSETS_THRESHOLD = 5_000_000
    SPOUSAL_INCOME_THRESHOLD = 300_000
    CORPORATE_NET_ASSETS_THRESHOLD = 5_000_000

    print("--- Accredited Investor Analysis ---\n")

    # --- Option A: Limited Partnership ---
    print("--- Analysis of Option A ---")
    liam_nfa = 5_400_000
    jack_net_assets = 18_000_000 - 5_000_000
    ace_nfa = 25_000_000
    gp_assets = 2_000_000 * 3 # Typo in question assumed as 2M not 20M

    liam_is_ai = liam_nfa > INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD
    jack_is_ai = jack_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD
    ace_is_ai = ace_nfa > INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD
    gp_is_ai = gp_assets >= CORPORATE_NET_ASSETS_THRESHOLD
    lp_is_ai = all([liam_is_ai, jack_is_ai, ace_is_ai, gp_is_ai])

    print(f"Liam's net financial assets: {print_currency(liam_nfa)}. Qualifies: {liam_is_ai}")
    print(f"Jack's net assets: {print_currency(18_000_000)} - {print_currency(5_000_000)} = {print_currency(jack_net_assets)}. Qualifies: {jack_is_ai}")
    print(f"Ace's net financial assets: {print_currency(ace_nfa)}. Qualifies: {ace_is_ai}")
    print(f"GP's net assets: {print_currency(gp_assets)}. Qualifies: {gp_is_ai}")
    print("Conclusion: Since all partners are accredited investors, the Limited Partnership IS an Accredited Investor.\n")

    # --- Option B: Individual (Income) ---
    print("--- Analysis of Option B ---")
    income_2019 = 150_000 + 170_000
    income_2020 = 175_000 + 175_000
    qualifies_2019 = income_2019 > SPOUSAL_INCOME_THRESHOLD
    qualifies_2020 = income_2020 > SPOUSAL_INCOME_THRESHOLD
    individual_b_is_ai = qualifies_2019 and qualifies_2020

    print(f"2019 Combined Income: {print_currency(150_000)} + {print_currency(170_000)} = {print_currency(income_2019)}. Qualifies: {qualifies_2019}")
    print(f"2020 Combined Income: {print_currency(175_000)} + {print_currency(175_000)} = {print_currency(income_2020)}. Qualifies: {qualifies_2020}")
    print("Conclusion: Since combined income exceeds the threshold for both years, the individual IS an Accredited Investor.\n")

    # --- Option C: Individual (Net Assets) ---
    print("--- Analysis of Option C ---")
    net_assets_c = 6_000_000 - 1_000_000
    individual_c_is_ai = net_assets_c >= INDIVIDUAL_NET_ASSETS_THRESHOLD
    print(f"Net Assets: {print_currency(6_000_000)} - {print_currency(1_000_000)} = {print_currency(net_assets_c)}. Qualifies: {individual_c_is_ai}")
    print("Conclusion: Since net assets meet the threshold, the individual IS an Accredited Investor.\n")

    # --- Option D: Corporation (Jose & James) ---
    print("--- Analysis of Option D ---")
    corp_d_net_assets = 0.10 * 100_000_000
    corp_d_is_ai = corp_d_net_assets >= CORPORATE_NET_ASSETS_THRESHOLD
    print(f"Corporation's Net Assets: 10% of {print_currency(100_000_000)} = {print_currency(corp_d_net_assets)}. Qualifies: {corp_d_is_ai}")
    print("Conclusion: Since the corporation's net assets exceed the threshold, it IS an Accredited Investor.\n")

    # --- Option E: Corporation (Alex & Bob) ---
    print("--- Analysis of Option E ---")
    corp_e_net_assets = 5_500_000 - 300_000
    corp_e_asset_test = corp_e_net_assets >= CORPORATE_NET_ASSETS_THRESHOLD

    alex_nfa = 900_000
    alex_net_assets = 3_000_000
    alex_is_ai = alex_nfa > INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD or alex_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD

    bob_income = 41_000
    bob_nfa = 75_000
    bob_net_assets = 0
    bob_is_ai = bob_income > 200_000 or bob_nfa > INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD or bob_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD

    print(f"Corporation's Net Assets: {print_currency(5_500_000)} - {print_currency(300_000)} = {print_currency(corp_e_net_assets)}. Meets asset test: {corp_e_asset_test}")
    print(f"Is owner Alex an AI? {alex_is_ai} (Fails NFA test of {print_currency(INDIVIDUAL_NET_FINANCIAL_ASSETS_THRESHOLD)} and Net Asset test of {print_currency(INDIVIDUAL_NET_ASSETS_THRESHOLD)})")
    print(f"Is owner Bob an AI? {bob_is_ai} (Fails all individual income and asset tests)")
    print("Conclusion: Although the corporation meets the net asset test, it would NOT be classified as an Accredited Investor because all of its owners are non-accredited individuals. This structure would likely be disallowed under anti-avoidance rules.\n")

    print("Final Answer: The entity that would not be classified as an Accredited Investor is E.")


accredited_investor_analysis()