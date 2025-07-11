def analyze_accredited_investors():
    """
    Analyzes several scenarios to determine which would not be classified as an Accredited Investor in Ontario.
    The analysis is based on National Instrument 45-106 as of January 2021.
    """

    # --- Accredited Investor Thresholds ---
    financial_assets_threshold = 1000000
    net_assets_threshold_individual = 5000000
    income_threshold_joint = 300000
    net_assets_threshold_entity = 5000000

    print("Analyzing each option based on Accredited Investor rules...\n")

    # --- Option A: Limited Partnership ---
    print("--- Option A: Limited Partnership ---")
    liam_net_financial_assets = 5400000
    jack_net_assets = 18000000 - 5000000
    ace_net_financial_assets = 25000000
    gp_inc_assets = 3 * 2000000 # Correcting typo in question 2,000,0000 -> 2,000,000
    
    liam_is_ai = liam_net_financial_assets > financial_assets_threshold
    jack_is_ai = jack_net_assets >= net_assets_threshold_individual
    ace_is_ai = ace_net_financial_assets > financial_assets_threshold
    gp_inc_is_ai = gp_inc_assets >= net_assets_threshold_entity
    
    print(f"Liam has ${liam_net_financial_assets:,.2f} in net financial assets, which is > ${financial_assets_threshold:,.2f}. Liam is an AI.")
    print(f"Jack has net assets of ${18000000:,.2f} - ${5000000:,.2f} = ${jack_net_assets:,.2f}, which is >= ${net_assets_threshold_individual:,.2f}. Jack is an AI.")
    print(f"Ace has ${ace_net_financial_assets:,.2f} in net financial assets, which is > ${financial_assets_threshold:,.2f}. Ace is an AI.")
    print(f"The GP, 2538901 Ontario Inc., has net assets of 3 * ${2000000:,.2f} = ${gp_inc_assets:,.2f}, which is >= ${net_assets_threshold_entity:,.2f}. The GP is an AI.")
    
    all_owners_are_ai_A = liam_is_ai and jack_is_ai and ace_is_ai and gp_inc_is_ai
    print("\nResult A: Since all limited and general partners are Accredited Investors, the Limited Partnership itself qualifies as an Accredited Investor.")
    print("-" * 40)

    # --- Option B: Individual (Joint Income Test) ---
    print("--- Option B: Individual (Joint Income Test) ---")
    joint_income_2019 = 150000 + 170000
    joint_income_2020 = 175000 + 175000
    
    print(f"2019 Joint Income: ${150000:,.2f} + ${170000:,.2f} = ${joint_income_2019:,.2f}")
    print(f"2020 Joint Income: ${175000:,.2f} + ${175000:,.2f} = ${joint_income_2020:,.2f}")
    
    is_ai_B = joint_income_2019 > income_threshold_joint and joint_income_2020 > income_threshold_joint
    print(f"\nResult B: Since joint income exceeded ${income_threshold_joint:,.2f} in both years, the individual is an Accredited Investor.")
    print("-" * 40)
    
    # --- Option C: Individual (Joint Net Asset Test) ---
    print("--- Option C: Individual (Joint Net Asset Test) ---")
    joint_net_assets = 6000000 - 1000000
    
    print(f"Joint Net Assets: ${6000000:,.2f} - ${1000000:,.2f} = ${joint_net_assets:,.2f}")
    
    is_ai_C = joint_net_assets >= net_assets_threshold_individual
    print(f"\nResult C: Since joint net assets of ${joint_net_assets:,.2f} meet the threshold of ${net_assets_threshold_individual:,.2f}, the individual is an Accredited Investor.")
    print("-" * 40)

    # --- Option D: Corporation (Jose and James) ---
    print("--- Option D: Corporation (Jose and James) ---")
    james_joint_income_2019 = 75000 + 210000
    james_joint_income_2020 = 75000 + 210000
    print(f"Checking owners: Jose has net financial assets of $100,000,000.00, making him an AI.")
    print(f"James's joint income in 2019 was ${75000:,.2f} + ${210000:,.2f} = ${james_joint_income_2019:,.2f}, which is below the ${income_threshold_joint:,.2f} threshold. James is not an AI.")
    print("Since not all owners are AIs, we test the corporation's assets.")
    
    corp_D_net_assets = 0.10 * 100000000
    print(f"Corporation's Net Assets: 10% of ${100000000:,.2f} = ${corp_D_net_assets:,.2f}")
    
    is_ai_D = corp_D_net_assets >= net_assets_threshold_entity
    print(f"\nResult D: Since the corporation's net assets of ${corp_D_net_assets:,.2f} exceed the ${net_assets_threshold_entity:,.2f} threshold, the corporation is an Accredited Investor.")
    print("-" * 40)

    # --- Option E: Corporation (Alex and Bob) ---
    print("--- Option E: Corporation (Alex and Bob) ---")
    alex_net_financial_assets = 900000
    alex_net_assets = 3000000
    print(f"Checking owners: Alex's net financial assets (${alex_net_financial_assets:,.2f}) are below ${financial_assets_threshold:,.2f} and his net assets (${alex_net_assets:,.2f}) are below ${net_assets_threshold_individual:,.2f}. Alex is not an AI.")
    print(f"Bob's income and assets do not meet any AI test thresholds. Bob is not an AI.")
    print("Since not all owners are AIs, we test the corporation's assets.")
    
    corp_E_net_assets = 5500000 - 300000
    print(f"Corporation's Net Assets: ${5500000:,.2f} - ${300000:,.2f} = ${corp_E_net_assets:,.2f}")
    
    # Mathematical result
    is_ai_E_math = corp_E_net_assets >= net_assets_threshold_entity
    print(f"Numerically, the corporation's net assets of ${corp_E_net_assets:,.2f} exceed the ${net_assets_threshold_entity:,.2f} threshold.")
    
    print("\nResult E: Although the corporation meets the numerical net asset test, securities regulators may not classify it as an AI.")
    print("This is because the entity is owned entirely by non-accredited individuals and appears to exist solely to pool their funds for investment, which can be viewed as a circumvention of prospectus requirements.")
    print("-" * 40)

    # --- Final Conclusion ---
    print("\nConclusion:")
    print("Options A, B, C, and D clearly qualify as Accredited Investors under standard tests.")
    print("Option E is the only one that is likely to NOT be classified as an Accredited Investor due to regulatory scrutiny of its structure, where non-accredited individuals have created an entity that numerically meets the test.")

analyze_accredited_investors()
<<<E>>>