def solve_accredited_investor_problem():
    """
    Analyzes five scenarios to determine which would not be classified as an Accredited Investor
    under Ontario securities laws as of January 2021.
    """

    # Accredited Investor Thresholds
    IND_FINANCIAL_ASSETS_THRESHOLD = 1000000
    IND_NET_ASSETS_THRESHOLD = 5000000
    IND_INCOME_THRESHOLD = 200000
    JOINT_INCOME_THRESHOLD = 300000
    CORP_NET_ASSETS_THRESHOLD = 5000000

    # --- Option A Analysis ---
    print("--- Analysis of Option A (Limited Partnership) ---")
    liam_nfa = 5400000
    jack_assets = 18000000
    jack_liabilities = 5000000
    ace_nfa = 25000000
    gp_assets_gifted = 2000000 * 3 # Gifted by Liam, Jack, and Ace to the GP

    liam_is_ai = liam_nfa > IND_FINANCIAL_ASSETS_THRESHOLD
    jack_net_assets = jack_assets - jack_liabilities
    jack_is_ai = jack_net_assets >= IND_NET_ASSETS_THRESHOLD
    ace_is_ai = ace_nfa > IND_FINANCIAL_ASSETS_THRESHOLD
    
    # The General Partner (GP) is an AI because its owners (Liam, Jack, Ace) are AIs,
    # and it also has net assets > $5M.
    gp_is_ai = (liam_is_ai and jack_is_ai and ace_is_ai) or (gp_assets_gifted >= CORP_NET_ASSETS_THRESHOLD)

    # The Limited Partnership (LP) is an AI if all its owners (LPs and GP) are AIs.
    lp_is_ai = liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai
    
    print(f"Liam's net financial assets: ${liam_nfa:,.2f}. Is Liam an AI? {liam_is_ai}")
    print(f"Jack's net assets: ${jack_assets:,.2f} - ${jack_liabilities:,.2f} = ${jack_net_assets:,.2f}. Is Jack an AI? {jack_is_ai}")
    print(f"Ace's net financial assets: ${ace_nfa:,.2f}. Is Ace an AI? {ace_is_ai}")
    print(f"The General Partner has owners who are all AIs. Is the GP an AI? {gp_is_ai}")
    print("Conclusion: All partners (limited and general) are Accredited Investors.")
    print("Therefore, the Limited Partnership in Option A is an Accredited Investor.\n")

    # --- Option B Analysis ---
    print("--- Analysis of Option B (Individual - Joint Income) ---")
    individual_income_2019 = 150000
    spouse_income_2019 = 170000
    individual_income_2020 = 175000
    spouse_income_2020 = 175000

    joint_income_2019 = individual_income_2019 + spouse_income_2019
    joint_income_2020 = individual_income_2020 + spouse_income_2020
    
    passes_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
    passes_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD
    individual_b_is_ai = passes_2019 and passes_2020

    print(f"Joint income in 2019: ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}")
    print(f"Joint income in 2020: ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}")
    print(f"Both years exceed the ${JOINT_INCOME_THRESHOLD:,.2f} threshold.")
    print("Conclusion: The individual in Option B is an Accredited Investor.\n")

    # --- Option C Analysis ---
    print("--- Analysis of Option C (Individual - Joint Net Assets) ---")
    joint_assets = 6000000
    joint_liabilities = 1000000
    
    joint_net_assets = joint_assets - joint_liabilities
    individual_c_is_ai = joint_net_assets >= IND_NET_ASSETS_THRESHOLD
    
    print(f"Joint net assets: ${joint_assets:,.2f} - ${joint_liabilities:,.2f} = ${joint_net_assets:,.2f}")
    print(f"The total meets the 'at least ${IND_NET_ASSETS_THRESHOLD:,.2f}' threshold.")
    print("Conclusion: The individual in Option C is an Accredited Investor.\n")

    # --- Option D Analysis ---
    print("--- Analysis of Option D (Corporation) ---")
    jose_nfa = 100000000
    transfer_pct = 0.10
    corp_d_assets = jose_nfa * transfer_pct
    
    corp_d_is_ai = corp_d_assets >= CORP_NET_ASSETS_THRESHOLD
    
    print(f"Corporation's net assets from transfer: {transfer_pct:.0%} of ${jose_nfa:,.2f} = ${corp_d_assets:,.2f}")
    print(f"The corporation's net assets exceed the ${CORP_NET_ASSETS_THRESHOLD:,.2f} threshold.")
    print("Conclusion: The corporation in Option D is an Accredited Investor on its own merit.\n")

    # --- Option E Analysis ---
    print("--- Analysis of Option E (Corporation) ---")
    corp_e_assets = 5500000
    corp_e_liabilities = 300000
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    
    alex_nfa = 900000
    alex_na = 3000000
    bob_income = 41000
    bob_nfa = 75000
    bob_na = 0
    
    # Check owners' status
    alex_is_ai = (alex_nfa > IND_FINANCIAL_ASSETS_THRESHOLD) or (alex_na >= IND_NET_ASSETS_THRESHOLD)
    bob_is_ai = (bob_income > IND_INCOME_THRESHOLD) or (bob_nfa > IND_FINANCIAL_ASSETS_THRESHOLD) or (bob_na >= IND_NET_ASSETS_THRESHOLD)
    
    # Check corporation's status
    corp_e_passes_net_asset_test = corp_e_net_assets >= CORP_NET_ASSETS_THRESHOLD
    corp_e_passes_owner_test = alex_is_ai and bob_is_ai

    print(f"Corporation's net assets: ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}")
    print(f"Does corporation pass the ${CORP_NET_ASSETS_THRESHOLD:,.2f} net asset test? {corp_e_passes_net_asset_test}")
    
    print("\nNow checking the status of the owners:")
    print(f"Alex's net financial assets (${alex_nfa:,.2f}) and net assets (${alex_na:,.2f}) are below the thresholds. Is Alex an AI? {alex_is_ai}")
    print(f"Bob's income, financial assets, and net assets are all below the thresholds. Is Bob an AI? {bob_is_ai}")
    print(f"Are all owners accredited investors? {corp_e_passes_owner_test}")
    
    print("\nConclusion: The corporation in Option E would not be classified as an Accredited Investor.")
    print("Although its net assets of $5,200,000.00 technically exceed the $5,000,000.00 threshold, it fails the 'all owners are AIs' test because neither Alex nor Bob are accredited investors. Securities regulations contain anti-avoidance provisions to prevent the use of entities by non-accredited individuals to circumvent investor protection laws. A corporation entirely owned and controlled by non-AIs, which only narrowly meets the asset test, is the most likely candidate to be disqualified upon review.")

    print("\n-------------------------------------------------")
    print("Final Answer: The option that would NOT be classified as an Accredited Investor is E.")
    print("-------------------------------------------------")

solve_accredited_investor_problem()
<<<E>>>