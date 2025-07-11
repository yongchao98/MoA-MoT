def solve():
    """
    Analyzes each option to determine if it qualifies as an Accredited Investor (AI)
    under Ontario securities laws as of January 2021.
    """
    # Define AI thresholds
    IND_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
    IND_NET_ASSETS_THRESHOLD = 5_000_000
    JOINT_INCOME_THRESHOLD = 300_000
    CORP_NET_ASSETS_THRESHOLD = 5_000_000

    print("--- Analysis of Accredited Investor Status ---\n")

    # --- Option A: Limited Partnership ---
    print("--- Analyzing Option A: Limited Partnership ---")
    liam_net_financial_assets = 5_400_000
    jack_net_assets = 18_000_000 - 5_000_000
    ace_net_financial_assets = 25_000_000
    gp_net_assets = 2_000_000 * 3
    
    liam_is_ai = liam_net_financial_assets > IND_FINANCIAL_ASSETS_THRESHOLD
    jack_is_ai = jack_net_assets >= IND_NET_ASSETS_THRESHOLD
    ace_is_ai = ace_net_financial_assets > IND_FINANCIAL_ASSETS_THRESHOLD
    gp_is_ai = gp_net_assets >= CORP_NET_ASSETS_THRESHOLD
    lp_is_ai = liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai

    print(f"Partner Liam's net financial assets of ${liam_net_financial_assets:,.2f} are > ${IND_FINANCIAL_ASSETS_THRESHOLD:,.2f}. Liam is an AI: {liam_is_ai}")
    print(f"Partner Jack's net assets are ${18_000_000:,.2f} - ${5_000_000:,.2f} = ${jack_net_assets:,.2f}, which is >= ${IND_NET_ASSETS_THRESHOLD:,.2f}. Jack is an AI: {jack_is_ai}")
    print(f"Partner Ace's net financial assets of ${ace_net_financial_assets:,.2f} are > ${IND_FINANCIAL_ASSETS_THRESHOLD:,.2f}. Ace is an AI: {ace_is_ai}")
    print(f"The General Partner corporation has net assets of 3 * ${2_000_000:,.2f} = ${gp_net_assets:,.2f}, which is >= ${CORP_NET_ASSETS_THRESHOLD:,.2f}. The GP is an AI: {gp_is_ai}")
    print("Conclusion: Since all partners are accredited investors, the limited partnership IS an Accredited Investor.\n")

    # --- Option B: Individual (Joint Income) ---
    print("--- Analyzing Option B: Individual (Joint Income Test) ---")
    joint_income_2019 = 150_000 + 170_000
    joint_income_2020 = 175_000 + 175_000
    individual_b_is_ai = joint_income_2019 > JOINT_INCOME_THRESHOLD and joint_income_2020 > JOINT_INCOME_THRESHOLD
    
    print(f"Couple's joint income in 2019 was ${150_000:,.2f} + ${170_000:,.2f} = ${joint_income_2019:,.2f}.")
    print(f"Couple's joint income in 2020 was ${175_000:,.2f} + ${175_000:,.2f} = ${joint_income_2020:,.2f}.")
    print(f"Both values are > the ${JOINT_INCOME_THRESHOLD:,.2f} threshold.")
    print("Conclusion: The individual IS an Accredited Investor.\n")
    
    # --- Option C: Individual (Joint Net Assets) ---
    print("--- Analyzing Option C: Individual (Joint Net Asset Test) ---")
    joint_net_assets = 6_000_000 - 1_000_000
    individual_c_is_ai = joint_net_assets >= IND_NET_ASSETS_THRESHOLD

    print(f"Couple's joint net assets are ${6_000_000:,.2f} - ${1_000_000:,.2f} = ${joint_net_assets:,.2f}.")
    print(f"This is >= the ${IND_NET_ASSETS_THRESHOLD:,.2f} threshold.")
    print("Conclusion: The individual IS an Accredited Investor.\n")

    # --- Option D: Corporation ---
    print("--- Analyzing Option D: Corporation ---")
    corp_d_assets = 0.10 * 100_000_000
    corp_d_is_ai = corp_d_assets >= CORP_NET_ASSETS_THRESHOLD
    james_joint_income = 75_000 + 210_000
    james_is_ai = james_joint_income > JOINT_INCOME_THRESHOLD

    print(f"The corporation's net assets are 10% of ${100_000_000:,.2f} = ${corp_d_assets:,.2f}.")
    print(f"This is > the ${CORP_NET_ASSETS_THRESHOLD:,.2f} threshold, so the corporation qualifies on its own merit.")
    print(f"(Note: Shareholder James' joint income of ${75_000:,.2f} + ${210_000:,.2f} = ${james_joint_income:,.2f} is below the ${JOINT_INCOME_THRESHOLD:,.2f} threshold, so he is not an AI. However, this does not disqualify the corporation since it qualifies on its own.)")
    print("Conclusion: The corporation IS an Accredited Investor.\n")
    
    # --- Option E: Corporation ---
    print("--- Analyzing Option E: Corporation ---")
    corp_e_net_assets_min = 5_500_000 - 300_000
    alex_is_ai = (900_000 > IND_FINANCIAL_ASSETS_THRESHOLD) or (3_000_000 >= IND_NET_ASSETS_THRESHOLD)
    bob_is_ai = (75_000 > IND_FINANCIAL_ASSETS_THRESHOLD) or (0 >= IND_NET_ASSETS_THRESHOLD)
    all_owners_ai = alex_is_ai and bob_is_ai

    print(f"The corporation's minimum net assets are ${5_500_000:,.2f} - ${300_000:,.2f} = ${corp_e_net_assets_min:,.2f}, which appears to meet the ${CORP_NET_ASSETS_THRESHOLD:,.2f} threshold.")
    print("However, let's check the owners:")
    print(f"Shareholder Alex has ${900_000:,.2f} in net financial assets and ${3_000_000:,.2f} in net assets. Both are below the AI thresholds. Alex is an AI: {alex_is_ai}.")
    print(f"Shareholder Bob has negligible assets and income. Bob is an AI: {bob_is_ai}.")
    print(f"Are all owners accredited investors? {all_owners_ai}")
    print("Conclusion: The corporation is owned entirely by non-accredited investors. Regulatory anti-avoidance rules prevent using such an entity to pool funds from non-AIs to qualify. Therefore, it would NOT be classified as an Accredited Investor.")
    
solve()
print("<<<E>>>")