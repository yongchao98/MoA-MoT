def analyze_accredited_investors():
    """
    Analyzes each option based on Ontario's Accredited Investor definitions
    as of January 2021 to determine which would not qualify.
    """

    print("Analyzing each option to determine its Accredited Investor status.\n")

    # --- Option A: Limited Partnership ---
    print("--- Analysis of Option A ---")
    liam_nfa = 5400000.00
    jack_net_assets = 18000000.00 - 5000000.00
    ace_nfa = 25000000.00
    gp_assets = 2000000.00 * 3

    liam_is_ai = liam_nfa > 1000000
    jack_is_ai = jack_net_assets >= 5000000
    ace_is_ai = ace_nfa > 1000000
    gp_is_ai = gp_assets > 5000000
    
    print(f"Liam's net financial assets of ${liam_nfa:,.2f} meet the >$1,000,000 test: {liam_is_ai}")
    print(f"Jack's net assets of ${18000000:,.2f} - ${5000000:,.2f} = ${jack_net_assets:,.2f} meet the >=$5,000,000 test: {jack_is_ai}")
    print(f"Ace's net financial assets of ${ace_nfa:,.2f} meet the >$1,000,000 test: {ace_is_ai}")
    print(f"The General Partner's assets of 3 * ${2000000:,.2f} = ${gp_assets:,.2f} meet the company >$5,000,000 NFA test: {gp_is_ai}")
    print("Conclusion: All partners are accredited investors. The limited partnership qualifies under the 'look-through' provision.")
    print("Status: IS an Accredited Investor.\n")


    # --- Option B: Individual (Joint Income) ---
    print("--- Analysis of Option B ---")
    joint_income_2019 = 150000.00 + 170000.00
    joint_income_2020 = 175000.00 + 175000.00
    is_ai = joint_income_2019 > 300000 and joint_income_2020 > 300000

    print(f"Couple's joint income in 2019: ${150000:,.2f} + ${170000:,.2f} = ${joint_income_2019:,.2f}")
    print(f"Couple's joint income in 2020: ${175000:,.2f} + ${175000:,.2f} = ${joint_income_2020:,.2f}")
    print(f"Both years' income exceed the $300,000 threshold: {is_ai}")
    print("Conclusion: The individual qualifies under the joint income test.")
    print("Status: IS an Accredited Investor.\n")


    # --- Option C: Individual (Joint Net Assets) ---
    print("--- Analysis of Option C ---")
    net_assets = 6000000.00 - 1000000.00
    is_ai = net_assets >= 5000000

    print(f"Couple's joint net assets: ${6000000:,.2f} - ${1000000:,.2f} = ${net_assets:,.2f}")
    print(f"Net assets meet the 'at least $5,000,000' test: {is_ai}")
    print("Conclusion: The individual qualifies under the joint net asset test.")
    print("Status: IS an Accredited Investor.\n")

    # --- Option D: Corporation (Funded by an AI) ---
    print("--- Analysis of Option D ---")
    corp_nfa = 100000000.00 * 0.10
    is_ai = corp_nfa > 5000000

    print(f"Corporation's net financial assets: ${100000000:,.2f} * 0.10 = ${corp_nfa:,.2f}")
    print(f"Corporation's net financial assets exceed the $5,000,000 test: {is_ai}")
    print("Conclusion: The corporation qualifies on its own, regardless of the status of its minority shareholder.")
    print("Status: IS an Accredited Investor.\n")

    # --- Option E: Corporation (Owned by non-AIs) ---
    print("--- Analysis of Option E ---")
    # First, check shareholders
    alex_is_ai = (900000.00 > 1000000) or (3000000.00 >= 5000000)
    bob_is_ai = (41000.00 > 200000) or (75000.00 > 1000000)
    print(f"Is shareholder Alex an accredited investor? {alex_is_ai}")
    print(f"Is shareholder Bob an accredited investor? {bob_is_ai}")
    print("Since not all shareholders are accredited investors, the corporation does not qualify via the 'look-through' test.")

    # Second, check the corporation itself
    corp_nfa = 5500000.00 - 300000.00
    is_ai_technically = corp_nfa > 5000000
    
    print(f"Corporation's net financial assets: ${5500000:,.2f} - ${300000:,.2f} = ${corp_nfa:,.2f}")
    print(f"Corporation's assets technically exceed the $5,000,000 test: {is_ai_technically}")
    
    print("Conclusion: Although technically passing the financial asset test, the corporation is entirely owned by non-accredited investors. This structure is likely to be considered a 'conduit' to pool funds and circumvent securities laws. Regulators would likely rule that it is not a legitimate accredited investor.")
    print("Status: Would NOT be classified as an Accredited Investor.\n")


analyze_accredited_investors()