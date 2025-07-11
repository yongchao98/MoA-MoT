def analyze_accredited_investors():
    """
    Analyzes several scenarios to determine which entity would not be classified
    as an Accredited Investor under Ontario securities laws as of Jan 2021.
    """
    
    # Define Accredited Investor (AI) thresholds
    INDIVIDUAL_INCOME_THRESHOLD = 200000
    JOINT_INCOME_THRESHOLD = 300000
    FINANCIAL_ASSETS_THRESHOLD = 1000000
    NET_ASSETS_THRESHOLD = 5000000
    CORPORATE_NET_ASSETS_THRESHOLD = 5000000

    print("--- Analyzing Accredited Investor Scenarios ---\n")

    # --- Option A: The Limited Partnership ---
    print("--- Option A: Limited Partnership Analysis ---")
    # Rule: A partnership is an AI if all of its owners are AIs.
    liam_net_financial_assets = 5400000.00
    jack_net_assets = 18000000.00 - 5000000.00
    ace_net_financial_assets = 25000000.00
    gp_corp_net_assets = 2000000.00 * 3 # Gifted amount

    is_liam_ai = liam_net_financial_assets > FINANCIAL_ASSETS_THRESHOLD
    print(f"Liam's net financial assets: ${liam_net_financial_assets:,.2f}. Test: > ${FINANCIAL_ASSETS_THRESHOLD:,.2f}. Liam is an AI: {is_liam_ai}")

    is_jack_ai = jack_net_assets >= NET_ASSETS_THRESHOLD
    print(f"Jack's net assets: ${18000000:,.2f} - ${5000000:,.2f} = ${jack_net_assets:,.2f}. Test: >= ${NET_ASSETS_THRESHOLD:,.2f}. Jack is an AI: {is_jack_ai}")

    is_ace_ai = ace_net_financial_assets > FINANCIAL_ASSETS_THRESHOLD
    print(f"Ace's net financial assets: ${ace_net_financial_assets:,.2f}. Test: > ${FINANCIAL_ASSETS_THRESHOLD:,.2f}. Ace is an AI: {is_ace_ai}")

    is_gp_corp_ai = gp_corp_net_assets >= CORPORATE_NET_ASSETS_THRESHOLD
    print(f"GP Corp's net assets: 3 * ${2000000:,.2f} = ${gp_corp_net_assets:,.2f}. Test: >= ${CORPORATE_NET_ASSETS_THRESHOLD:,.2f}. GP Corp is an AI: {is_gp_corp_ai}")

    is_lp_ai = all([is_liam_ai, is_jack_ai, is_ace_ai, is_gp_corp_ai])
    print(f"Conclusion: All owners are accredited investors. Therefore, the Limited Partnership IS an accredited investor.\n")


    # --- Option B: The Individual (Income) ---
    print("--- Option B: Individual Income Analysis ---")
    # Rule: Joint income > $300k for the last 2 years.
    joint_income_2019 = 150000.00 + 170000.00
    joint_income_2020 = 175000.00 + 175000.00
    
    meets_2019_income = joint_income_2019 > JOINT_INCOME_THRESHOLD
    print(f"2019 Joint Income: ${150000:,.2f} + ${170000:,.2f} = ${joint_income_2019:,.2f}. Test: > ${JOINT_INCOME_THRESHOLD:,.2f}. Met: {meets_2019_income}")
    
    meets_2020_income = joint_income_2020 > JOINT_INCOME_THRESHOLD
    print(f"2020 Joint Income: ${175000:,.2f} + ${175000:,.2f} = ${joint_income_2020:,.2f}. Test: > ${JOINT_INCOME_THRESHOLD:,.2f}. Met: {meets_2020_income}")

    is_b_ai = meets_2019_income and meets_2020_income
    print(f"Conclusion: The joint income test is passed for both years. Therefore, the individual IS an accredited investor.\n")

    # --- Option C: The Individual (Net Assets) ---
    print("--- Option C: Individual Net Asset Analysis ---")
    # Rule: Joint net assets >= $5M.
    joint_net_assets = 6000000.00 - 1000000.00
    
    is_c_ai = joint_net_assets >= NET_ASSETS_THRESHOLD
    print(f"Joint Net Assets: ${6000000:,.2f} - ${1000000:,.2f} = ${joint_net_assets:,.2f}. Test: >= ${NET_ASSETS_THRESHOLD:,.2f}. Met: {is_c_ai}")
    print(f"Conclusion: The joint net asset test is passed. Therefore, the individual IS an accredited investor.\n")
    
    # --- Option D: The Corporation (Jose and James) ---
    print("--- Option D: Corporation (Jose/James) Analysis ---")
    # Rule: A corporation is an AI if it has net assets >= $5M.
    corp_d_net_assets = 100000000.00 * 0.10
    
    is_corp_d_ai = corp_d_net_assets >= CORPORATE_NET_ASSETS_THRESHOLD
    print(f"Corporation's Net Assets: 10% of ${100000000:,.2f} = ${corp_d_net_assets:,.2f}. Test: >= ${CORPORATE_NET_ASSETS_THRESHOLD:,.2f}. Met: {is_corp_d_ai}")

    # Check James just to confirm the 'look-through' rule wouldn't apply
    james_joint_income = 75000.00 + 210000.00
    is_james_ai = james_joint_income > JOINT_INCOME_THRESHOLD
    print(f"Note: Co-owner James's joint income is ${75000:,.2f} + ${210000:,.2f} = ${james_joint_income:,.2f}, which is less than the ${JOINT_INCOME_THRESHOLD:,.2f} threshold. He is not an AI.")
    print("Conclusion: Although not all owners are AIs, the corporation qualifies on its own with over $5M in net assets. Therefore, the corporation IS an accredited investor.\n")
    
    # --- Option E: The Corporation (Alex and Bob) ---
    print("--- Option E: Corporation (Alex/Bob) Analysis ---")
    # Rule: The $5M net asset test excludes 'investment funds'. This corporation's sole purpose seems to be investing pooled money, so it is likely an investment fund.
    # Therefore, it must qualify based on its owners, who must all be accredited investors.
    print("Logic: This corporation's assets (GIC, investment properties) and structure suggest it's an 'investment fund', which is excluded from the standard $5M net asset test for corporations.")
    print("For an investment fund structured this way to be an AI, all its owners must be AIs.")
    
    # Analyze owner Alex
    alex_net_financial_assets = 900000.00
    alex_net_assets = 3000000.00
    is_alex_ai = (alex_net_financial_assets > FINANCIAL_ASSETS_THRESHOLD) or (alex_net_assets >= NET_ASSETS_THRESHOLD)
    print(f"Owner Alex's net financial assets: ${alex_net_financial_assets:,.2f} (Fails > ${FINANCIAL_ASSETS_THRESHOLD:,.2f} test).")
    print(f"Owner Alex's total net assets: ${alex_net_assets:,.2f} (Fails >= ${NET_ASSETS_THRESHOLD:,.2f} test). Alex is an AI: {is_alex_ai}")
    
    # Analyze owner Bob
    bob_income_2020 = 41000.00
    bob_financial_assets = 75000.00
    is_bob_ai = (bob_income_2020 > INDIVIDUAL_INCOME_THRESHOLD) or (bob_financial_assets > FINANCIAL_ASSETS_THRESHOLD)
    print(f"Owner Bob's income: ${bob_income_2020:,.2f} (Fails > ${INDIVIDUAL_INCOME_THRESHOLD:,.2f} test).")
    print(f"Owner Bob's financial assets: ${bob_financial_assets:,.2f} (Fails > ${FINANCIAL_ASSETS_THRESHOLD:,.2f} test). Bob is an AI: {is_bob_ai}")

    is_corp_e_ai = is_alex_ai and is_bob_ai
    print(f"Conclusion: Neither Alex nor Bob are accredited investors. Because the corporation is an investment fund and its owners are not AIs, it cannot qualify. Therefore, this corporation would NOT be classified as an accredited investor.\n")

    print("\n--- Final Determination ---")
    print("Based on the analysis, the entity in Option E is the one that would not be classified as an accredited investor.")
    print("<<<E>>>")

if __name__ == '__main__':
    analyze_accredited_investors()