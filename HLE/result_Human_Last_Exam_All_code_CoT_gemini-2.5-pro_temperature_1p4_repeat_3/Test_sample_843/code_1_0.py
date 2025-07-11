def solve_accredited_investor_case():
    """
    Analyzes several scenarios to determine which entity would not be classified as an Accredited Investor
    in Ontario as of January 2021.
    """

    # --- Accredited Investor Thresholds ---
    INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
    INDIVIDUAL_NET_ASSETS_THRESHOLD = 5_000_000
    INDIVIDUAL_INCOME_THRESHOLD = 200_000
    JOINT_INCOME_THRESHOLD = 300_000
    CORPORATE_NET_ASSETS_THRESHOLD = 5_000_000

    print("Analyzing each option based on Accredited Investor (AI) definitions in Ontario...\n")

    # --- Option A: Limited Partnership ---
    print("--- Option A: Limited Partnership ---")
    liam_net_financial_assets = 5_400_000
    jack_total_assets = 18_000_000
    jack_total_liabilities = 5_000_000
    ace_net_financial_assets = 25_000_000
    gp_gift = 2_000_000 * 3

    jack_net_assets = jack_total_assets - jack_total_liabilities
    
    liam_is_ai = liam_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD
    jack_is_ai = jack_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD
    ace_is_ai = ace_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD
    gp_is_ai = gp_gift >= CORPORATE_NET_ASSETS_THRESHOLD

    print(f"Liam's net financial assets are ${liam_net_financial_assets:,.2f}. Is Liam an AI? {liam_is_ai}")
    print(f"Jack's net assets are ${jack_total_assets:,.2f} - ${jack_total_liabilities:,.2f} = ${jack_net_assets:,.2f}. Is Jack an AI? {jack_is_ai}")
    print(f"Ace's net financial assets are ${ace_net_financial_assets:,.2f}. Is Ace an AI? {ace_is_ai}")
    print(f"The General Partner's net assets are at least ${gp_gift:,.2f}. Is the GP an AI? {gp_is_ai}")
    
    partnership_is_ai = liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai
    print("Conclusion: Since all owners are AIs, the Limited Partnership qualifies as an AI under the look-through provision.")
    print(f"Result for A: IS an Accredited Investor.\n")


    # --- Option B: Individual (Joint Income) ---
    print("--- Option B: Individual (Joint Income) ---")
    individual_income_2019 = 150_000
    spouse_income_2019 = 170_000
    individual_income_2020 = 175_000
    spouse_income_2020 = 175_000
    
    joint_income_2019 = individual_income_2019 + spouse_income_2019
    joint_income_2020 = individual_income_2020 + spouse_income_2020
    
    meets_2019_test = joint_income_2019 > JOINT_INCOME_THRESHOLD
    meets_2020_test = joint_income_2020 > JOINT_INCOME_THRESHOLD
    
    print(f"Combined 2019 income: ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}. Exceeds ${JOINT_INCOME_THRESHOLD:,.2f}? {meets_2019_test}")
    print(f"Combined 2020 income: ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}. Exceeds ${JOINT_INCOME_THRESHOLD:,.2f}? {meets_2020_test}")
    
    individual_b_is_ai = meets_2019_test and meets_2020_test
    print("Conclusion: The individual qualifies as an AI based on the joint income test.")
    print(f"Result for B: IS an Accredited Investor.\n")
    

    # --- Option C: Individual (Joint Net Assets) ---
    print("--- Option C: Individual (Joint Net Assets) ---")
    total_assets_c = 6_000_000
    total_liabilities_c = 1_000_000
    net_assets_c = total_assets_c - total_liabilities_c
    
    individual_c_is_ai = net_assets_c >= INDIVIDUAL_NET_ASSETS_THRESHOLD
    
    print(f"Combined net assets: ${total_assets_c:,.2f} - ${total_liabilities_c:,.2f} = ${net_assets_c:,.2f}.")
    print(f"Is net assets >= ${INDIVIDUAL_NET_ASSETS_THRESHOLD:,.2f}? {individual_c_is_ai}")
    print("Conclusion: The individual qualifies as an AI based on the joint net asset test.")
    print(f"Result for C: IS an Accredited Investor.\n")

    # --- Option D: Corporation (Jose and James) ---
    print("--- Option D: Corporation (Jose and James) ---")
    jose_nfa = 100_000_000
    transfer_pct = 0.10
    corp_d_assets = jose_nfa * transfer_pct
    
    corp_d_is_ai = corp_d_assets >= CORPORATE_NET_ASSETS_THRESHOLD
    print("Test 1: Corporation's Net Assets")
    print(f"Assets transferred from Jose: {transfer_pct:.0%} of ${jose_nfa:,.2f} = ${corp_d_assets:,.2f}.")
    print(f"Are corporate net assets >= ${CORPORATE_NET_ASSETS_THRESHOLD:,.2f}? {corp_d_is_ai}")
    print("Conclusion: The corporation qualifies as an AI based on its own net assets.")
    print(f"Result for D: IS an Accredited Investor.\n")

    # --- Option E: Corporation (Alex and Bob) ---
    print("--- Option E: Corporation (Alex and Bob) ---")
    corp_e_financial_assets = 5_500_000
    corp_e_financial_liabilities = 300_000
    corp_e_net_financial_assets = corp_e_financial_assets - corp_e_financial_liabilities
    
    alex_nfa = 900_000
    alex_net_assets = 3_000_000
    bob_income = 41_000
    
    alex_is_ai = (alex_nfa > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD) or (alex_net_assets >= INDIVIDUAL_NET_ASSETS_THRESHOLD)
    bob_is_ai = False # Income and assets are far below thresholds
    
    look_through_test_passed = alex_is_ai and bob_is_ai
    
    print("Test 1: Corporation's Net Assets")
    print(f"The corporation's total net assets cannot be determined. The value of its 'investment properties' is unknown. The test for a corporation is based on TOTAL net assets, not just financial assets. Therefore, we cannot confirm it meets the ${CORPORATE_NET_ASSETS_THRESHOLD:,.2f} test.")
    
    print("\nTest 2: Look-through Provision (Are all owners AIs?)")
    print(f"Is Alex an AI? Alex's net financial assets (${alex_nfa:,.2f}) are below ${INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD:,.2f} and net assets (${alex_net_assets:,.2f}) are below ${INDIVIDUAL_NET_ASSETS_THRESHOLD:,.2f}. Status: {alex_is_ai}")
    print(f"Is Bob an AI? Bob's income and assets are below all thresholds. Status: {bob_is_ai}")
    print(f"Does the corporation pass the look-through test? {look_through_test_passed}")
    
    print("\nConclusion: The corporation fails the look-through test because neither Alex nor Bob is an AI. Furthermore, there is insufficient information to prove it meets the $5,000,000 net asset test.")
    print(f"Result for E: Would NOT be classified as an Accredited Investor.\n")

solve_accredited_investor_case()