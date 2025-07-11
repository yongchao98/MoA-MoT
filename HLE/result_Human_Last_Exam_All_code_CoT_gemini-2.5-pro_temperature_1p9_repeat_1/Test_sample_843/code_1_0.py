def solve():
    """
    Analyzes five different scenarios to determine which would not be classified
    as an Accredited Investor under Ontario securities laws as of January 2021.
    """
    
    # Define Accredited Investor Thresholds
    FINANCIAL_ASSETS_THRESHOLD_IND = 1000000
    NET_ASSETS_THRESHOLD_IND = 5000000
    JOINT_INCOME_THRESHOLD = 300000
    NET_ASSETS_THRESHOLD_ENTITY = 5000000

    print("--- Analysis of Accredited Investor Status ---")

    # --- Option A: Limited Partnership ---
    print("\n[Option A: Limited Partnership Analysis]")
    liam_nfa = 5400000.00
    is_liam_accredited = liam_nfa > FINANCIAL_ASSETS_THRESHOLD_IND
    print(f"Liam's net financial assets: ${liam_nfa:,.2f}")
    print(f"Is Liam's NFA > ${FINANCIAL_ASSETS_THRESHOLD_IND:,.2f}? {is_liam_accredited}. Liam is accredited.")

    jack_assets = 18000000.00
    jack_liabilities = 5000000.00
    jack_net_assets = jack_assets - jack_liabilities
    is_jack_accredited = jack_net_assets > NET_ASSETS_THRESHOLD_IND
    print(f"Jack's net assets: ${jack_assets:,.2f} - ${jack_liabilities:,.2f} = ${jack_net_assets:,.2f}")
    print(f"Is Jack's Net Assets > ${NET_ASSETS_THRESHOLD_IND:,.2f}? {is_jack_accredited}. Jack is accredited.")
    
    ace_nfa = 25000000.00
    is_ace_accredited = ace_nfa > FINANCIAL_ASSETS_THRESHOLD_IND
    print(f"Ace's net financial assets: ${ace_nfa:,.2f}")
    print(f"Is Ace's NFA > ${FINANCIAL_ASSETS_THRESHOLD_IND:,.2f}? {is_ace_accredited}. Ace is accredited.")
    print("The general partner is owned by Liam, Jack, and Ace, who are all accredited. Thus the GP is accredited.")
    print("Result for A: The Limited Partnership is ACCREDITED because all of its owners are accredited.")

    # --- Option B: Individual (Joint Income) ---
    print("\n[Option B: Individual - Joint Income Analysis]")
    ind_income_2019 = 150000.00
    spouse_income_2019 = 170000.00
    joint_income_2019 = ind_income_2019 + spouse_income_2019
    is_2019_met = joint_income_2019 > JOINT_INCOME_THRESHOLD
    print(f"2019 Joint Income: ${ind_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}")
    print(f"Is 2019 Joint Income > ${JOINT_INCOME_THRESHOLD:,.2f}? {is_2019_met}")
    
    ind_income_2020 = 175000.00
    spouse_income_2020 = 175000.00
    joint_income_2020 = ind_income_2020 + spouse_income_2020
    is_2020_met = joint_income_2020 > JOINT_INCOME_THRESHOLD
    print(f"2020 Joint Income: ${ind_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}")
    print(f"Is 2020 Joint Income > ${JOINT_INCOME_THRESHOLD:,.2f}? {is_2020_met}")
    print("Result for B: The individual is ACCREDITED based on the joint income test.")

    # --- Option C: Individual (Joint Net Assets) ---
    print("\n[Option C: Individual - Joint Net Assets Analysis]")
    total_assets = 6000000.00
    total_liabilities = 1000000.00
    net_assets = total_assets - total_liabilities
    # The rule requires net assets to EXCEED the threshold
    is_met = net_assets > NET_ASSETS_THRESHOLD_IND
    print(f"Joint Net Assets: ${total_assets:,.2f} - ${total_liabilities:,.2f} = ${net_assets:,.2f}")
    print(f"Is Joint Net Assets > ${NET_ASSETS_THRESHOLD_IND:,.2f}? {is_met}")
    print("Result for C: The individual is NOT ACCREDITED because their net assets do not exceed $5,000,000.")

    # --- Option D: Corporation (Jose and James) ---
    print("\n[Option D: Corporation Analysis]")
    # A corporation is accredited if it meets EITHER the entity net asset test OR the all-owners-accredited test.
    jose_nfa = 100000000.00
    corp_assets = 0.10 * jose_nfa
    corp_liabilities = 0
    corp_net_assets = corp_assets - corp_liabilities
    is_corp_accredited_by_assets = corp_net_assets >= NET_ASSETS_THRESHOLD_ENTITY
    print(f"Corporation's Assets: 10% of ${jose_nfa:,.2f} = ${corp_assets:,.2f}")
    print(f"Corporation's Net Assets = ${corp_net_assets:,.2f}")
    print(f"Is Corporation's Net Assets >= ${NET_ASSETS_THRESHOLD_ENTITY:,.2f}? {is_corp_accredited_by_assets}")
    print("Result for D: The Corporation is ACCREDITED because its net assets are at least $5,000,000.")
    
    # --- Option E: Corporation (Alex and Bob) ---
    print("\n[Option E: Corporation Analysis]")
    corp_e_assets = 5500000.00
    corp_e_liabilities = 300000.00
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    is_corp_e_accredited = corp_e_net_assets >= NET_ASSETS_THRESHOLD_ENTITY
    print(f"Corporation's Net Assets: ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}")
    print(f"Is Corporation's Net Assets >= ${NET_ASSETS_THRESHOLD_ENTITY:,.2f}? {is_corp_e_accredited}")
    print("Result for E: The Corporation is ACCREDITED because its net assets are at least $5,000,000.")

solve()
print("\n<<<C>>>")