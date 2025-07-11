def solve():
    """
    Analyzes each option to determine if it qualifies as an Accredited Investor (AI)
    under Ontario securities laws as of January 2021.
    """
    
    # Define Accredited Investor thresholds in CAD
    IND_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
    JOINT_INCOME_THRESHOLD = 300_000
    IND_NET_ASSETS_THRESHOLD = 5_000_000
    CORP_NET_ASSETS_THRESHOLD = 5_000_000

    print("--- Analysis of Each Option ---")

    # --- Option A: The Limited Partnership ---
    print("\n[A] Analyzing the Limited Partnership (LP):")
    print("An LP qualifies if all its partners are Accredited Investors (AIs).")
    liam_net_financial_assets = 5_400_000
    jack_net_assets = 18_000_000 - 5_000_000
    ace_net_financial_assets = 25_000_000
    # The General Partner received 3 gifts of 2,000,000
    gp_net_assets = 3 * 2_000_000

    liam_is_ai = liam_net_financial_assets > IND_FINANCIAL_ASSETS_THRESHOLD
    jack_is_ai = jack_net_assets >= IND_NET_ASSETS_THRESHOLD
    ace_is_ai = ace_net_financial_assets > IND_FINANCIAL_ASSETS_THRESHOLD
    gp_is_ai = gp_net_assets >= CORP_NET_ASSETS_THRESHOLD
    
    print(f"Liam's net financial assets are ${liam_net_financial_assets:,.2f}. Qualifies as AI: {liam_is_ai}")
    print(f"Jack's net assets are $18,000,000 - $5,000,000 = ${jack_net_assets:,.2f}. Qualifies as AI: {jack_is_ai}")
    print(f"Ace's net financial assets are ${ace_net_financial_assets:,.2f}. Qualifies as AI: {ace_is_ai}")
    print(f"The General Partner's net assets are ${gp_net_assets:,.2f}. Qualifies as AI: {gp_is_ai}")
    
    if liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai:
        print("Conclusion: All partners are AIs, therefore the Limited Partnership IS an Accredited Investor.")
    else:
        print("Conclusion: Not all partners are AIs, therefore the Limited Partnership IS NOT an Accredited Investor.")

    # --- Option B: The Individual (Joint Income) ---
    print("\n[B] Analyzing the individual based on joint income:")
    print("An individual qualifies if their joint income with a spouse exceeded $300,000 in each of the past two years.")
    individual_income_2019 = 150_000
    spouse_income_2019 = 170_000
    joint_income_2019 = individual_income_2019 + spouse_income_2019

    individual_income_2020 = 175_000
    spouse_income_2020 = 175_000
    joint_income_2020 = individual_income_2020 + spouse_income_2020
    
    income_test_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
    income_test_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD
    
    print(f"Joint income in 2019 was ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}. Met threshold: {income_test_2019}")
    print(f"Joint income in 2020 was ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}. Met threshold: {income_test_2020}")

    if income_test_2019 and income_test_2020:
        print("Conclusion: The joint income test is met for both years. The individual IS an Accredited Investor.")
    else:
        print("Conclusion: The joint income test was not met. The individual IS NOT an Accredited Investor.")

    # --- Option C: The Individual (Joint Net Assets) ---
    print("\n[C] Analyzing the individual based on joint net assets:")
    print("An individual qualifies if their joint net assets with a spouse are at least $5,000,000.")
    joint_assets = 6_000_000
    joint_liabilities = 1_000_000
    joint_net_assets = joint_assets - joint_liabilities

    net_asset_test = joint_net_assets >= IND_NET_ASSETS_THRESHOLD
    
    print(f"Joint net assets are ${joint_assets:,.2f} - ${joint_liabilities:,.2f} = ${joint_net_assets:,.2f}.")
    print(f"The assets meet the 'at least ${IND_NET_ASSETS_THRESHOLD:,.2f}' threshold: {net_asset_test}")
    if net_asset_test:
        print("Conclusion: The joint net asset test is met. The individual IS an Accredited Investor.")
    else:
        print("Conclusion: The joint net asset test is not met. The individual IS NOT an Accredited Investor.")

    # --- Option D: The Corporation (Jose and James) ---
    print("\n[D] Analyzing the corporation owned by Jose and James:")
    print("A corporation qualifies if its own net assets are at least $5,000,000.")
    jose_net_financial_assets = 100_000_000
    transfer_percentage = 0.10
    assets_transferred = jose_net_financial_assets * transfer_percentage
    
    corp_d_net_assets_test = assets_transferred >= CORP_NET_ASSETS_THRESHOLD

    print(f"The corporation's net assets are ${assets_transferred:,.2f} (from Jose's transfer).")
    print(f"The corporation's assets meet the '${CORP_NET_ASSETS_THRESHOLD:,.2f}' threshold: {corp_d_net_assets_test}")
    if corp_d_net_assets_test:
        print("Conclusion: The corporation's net asset test is met. The corporation IS an Accredited Investor.")
    else:
        print("Conclusion: The corporation's net asset test is not met.")

    # --- Option E: The Corporation (Alex and Bob) ---
    print("\n[E] Analyzing the corporation owned by Alex and Bob:")
    print("Checking if the corporation or its owners qualify.")
    
    # Check owners first
    alex_net_financial_assets = 900_000
    bob_income_2020 = 41_000
    alex_is_ai = alex_net_financial_assets > IND_FINANCIAL_ASSETS_THRESHOLD
    bob_is_ai = False # Fails income and asset tests
    print(f"Alex's net financial assets are ${alex_net_financial_assets:,.2f}. Qualifies as AI: {alex_is_ai}")
    print(f"Bob's reported income and assets do not meet any AI test. Qualifies as AI: {bob_is_ai}")
    print("Since not all owners are AIs, the look-through test fails. Now check the corporation's own net assets.")
    
    # Check corporation's net assets
    corp_e_assets = 5_500_000
    corp_e_liabilities = 300_000
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    
    corp_e_net_asset_test = corp_e_net_assets >= CORP_NET_ASSETS_THRESHOLD
    
    print(f"The corporation's stated net assets are ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}.")
    print(f"On paper, the assets meet the '${CORP_NET_ASSETS_THRESHOLD:,.2f}' threshold: {corp_e_net_asset_test}")
    
    print("\nFinal consideration for [E]:")
    print("Despite the calculation, securities regulations include anti-avoidance provisions.")
    print("A corporation created by non-accredited individuals (like Alex and Bob) solely to pool assets and meet the $5,000,000 threshold may not be recognized by regulators.")
    print("Given that the owners are not AIs and the company's main asset is passive (a GIC), it is likely to be viewed as such a vehicle.")
    print("Conclusion: The corporation would likely NOT be classified as an Accredited Investor.")


solve()
<<<E>>>