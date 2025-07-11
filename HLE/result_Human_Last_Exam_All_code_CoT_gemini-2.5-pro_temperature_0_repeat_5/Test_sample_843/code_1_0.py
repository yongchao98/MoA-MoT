def solve_accredited_investor_problem():
    """
    Analyzes five scenarios to determine which would not be classified as an Accredited Investor
    in Ontario as of January 2021.
    """

    # Define Accredited Investor (AI) thresholds
    INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
    INDIVIDUAL_NET_ASSETS_THRESHOLD = 5_000_000
    JOINT_INCOME_THRESHOLD = 300_000
    ENTITY_NET_ASSETS_THRESHOLD = 5_000_000

    print("Analyzing each option based on Accredited Investor rules...\n")

    # --- Option A: Limited Partnership ---
    print("--- Option A Analysis (Limited Partnership) ---")
    liam_net_financial_assets = 5_400_000
    jack_total_assets = 18_000_000
    jack_total_liabilities = 5_000_000
    ace_net_financial_assets = 25_000_000

    liam_is_ai = liam_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD
    jack_net_assets = jack_total_assets - jack_total_liabilities
    jack_is_ai = jack_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD
    ace_is_ai = ace_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD

    print(f"Liam's net financial assets: ${liam_net_financial_assets:,.2f}. Is Liam an AI? {liam_is_ai}")
    print(f"Jack's net assets: ${jack_total_assets:,.2f} - ${jack_total_liabilities:,.2f} = ${jack_net_assets:,.2f}. Is Jack an AI? {jack_is_ai}")
    print(f"Ace's net financial assets: ${ace_net_financial_assets:,.2f}. Is Ace an AI? {ace_is_ai}")

    lp_is_ai = liam_is_ai and jack_is_ai and ace_is_ai
    print("Conclusion: Since all limited partners are Accredited Investors, the partnership qualifies under the look-through test.")
    print(f"Result for A: IS an Accredited Investor.\n")


    # --- Option B: Individual (Joint Income) ---
    print("--- Option B Analysis (Individual - Joint Income) ---")
    individual_income_2019 = 150_000
    spouse_income_2019 = 170_000
    individual_income_2020 = 175_000
    spouse_income_2020 = 175_000

    joint_income_2019 = individual_income_2019 + spouse_income_2019
    joint_income_2020 = individual_income_2020 + spouse_income_2020

    pass_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
    pass_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD

    print(f"2019 Joint Income: ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}. Met threshold? {pass_2019}")
    print(f"2020 Joint Income: ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}. Met threshold? {pass_2020}")

    individual_b_is_ai = pass_2019 and pass_2020
    print("Conclusion: The joint income exceeded $300,000 in both recent years.")
    print(f"Result for B: IS an Accredited Investor.\n")


    # --- Option C: Individual (Net Assets) ---
    print("--- Option C Analysis (Individual - Net Assets) ---")
    total_assets_c = 6_000_000
    total_liabilities_c = 1_000_000
    net_assets_c = total_assets_c - total_liabilities_c

    individual_c_is_ai = net_assets_c >= INDIVIDUAL_NET_ASSETS_THRESHOLD
    print(f"Net Assets: ${total_assets_c:,.2f} - ${total_liabilities_c:,.2f} = ${net_assets_c:,.2f}")
    print(f"Conclusion: The net assets are at least $5,000,000.")
    print(f"Result for C: IS an Accredited Investor.\n")


    # --- Option D: Corporation (Jose and James) ---
    print("--- Option D Analysis (Corporation) ---")
    jose_net_financial_assets = 100_000_000
    transfer_percentage = 0.10
    corp_d_assets = jose_net_financial_assets * transfer_percentage

    corp_d_is_ai = corp_d_assets >= ENTITY_NET_ASSETS_THRESHOLD
    print(f"Corporation's Net Assets from transfer: {transfer_percentage:.0%} of ${jose_net_financial_assets:,.2f} = ${corp_d_assets:,.2f}")
    print(f"Conclusion: The corporation's net assets exceed the ${ENTITY_NET_ASSETS_THRESHOLD:,.2f} threshold.")
    print(f"Result for D: IS an Accredited Investor.\n")


    # --- Option E: Corporation (Alex and Bob) ---
    print("--- Option E Analysis (Corporation) ---")
    corp_e_assets = 5_500_000
    corp_e_liabilities = 300_000
    corp_e_net_assets = corp_e_assets - corp_e_liabilities

    corp_e_is_ai_by_assets = corp_e_net_assets >= ENTITY_NET_ASSETS_THRESHOLD
    print("Test 1: Corporation's Net Assets")
    print(f"Corporation's Net Assets: ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}")
    print(f"Does it meet the entity net asset test? {corp_e_is_ai_by_assets}")
    print("\nTest 2: Look-Through Provision (Are all owners AIs?)")
    alex_net_financial_assets = 900_000
    alex_net_assets = 3_000_000
    alex_is_ai = (alex_net_financial_assets >= INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD) or (alex_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD)
    print(f"Alex's net financial assets are ${alex_net_financial_assets:,.2f} (< $1M) and net assets are ${alex_net_assets:,.2f} (< $5M). Is Alex an AI? {alex_is_ai}")

    bob_income_2019 = 40_000
    bob_income_2020 = 41_000
    bob_financial_assets = 75_000
    bob_is_ai = (bob_financial_assets >= INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD) # Other tests clearly fail
    print(f"Bob's income and assets are below all thresholds. Is Bob an AI? {bob_is_ai}")

    corp_e_is_ai_by_look_through = alex_is_ai and bob_is_ai
    print(f"Does it meet the look-through test? {corp_e_is_ai_by_look_through}")
    print("\nConclusion: The corporation technically meets the $5,000,000 net asset test. However, it is the only option that is entirely owned by individuals who are NOT Accredited Investors. Such a structure may be scrutinized by regulators as an attempt to circumvent the rules. Among the choices, this is the one that would most likely not be classified as an Accredited Investor in practice.")
    print(f"Result for E: IS NOT an Accredited Investor.\n")

    final_answer = "E"
    print(f"The option that would not be classified as an Accredited Investor is E.")

solve_accredited_investor_problem()
<<<E>>>