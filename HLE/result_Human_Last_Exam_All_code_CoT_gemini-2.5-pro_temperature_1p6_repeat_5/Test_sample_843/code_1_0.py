def solve_accredited_investor_case():
    """
    Analyzes several scenarios to determine which entity would not be classified
    as an Accredited Investor in Ontario as of January 2021.
    """

    print("--- Analyzing Accredited Investor Status for Each Option ---\n")

    # --- Option A: Limited Partnership ---
    print("A: Checking Limited Partnership")
    liam_net_financial_assets = 5_400_000
    liam_is_ai = liam_net_financial_assets > 1_000_000
    print(f"  - Liam's net financial assets: ${liam_net_financial_assets:,.2f}. Meets the > $1,000,000 test: {liam_is_ai}")

    jack_assets = 18_000_000
    jack_liabilities = 5_000_000
    jack_net_assets = jack_assets - jack_liabilities
    jack_is_ai = jack_net_assets >= 5_000_000
    print(f"  - Jack's net assets: ${jack_assets:,.2f} - ${jack_liabilities:,.2f} = ${jack_net_assets:,.2f}. Meets the >= $5,000,000 test: {jack_is_ai}")

    ace_net_financial_assets = 25_000_000
    ace_is_ai = ace_net_financial_assets > 1_000_000
    print(f"  - Ace's net financial assets: ${ace_net_financial_assets:,.2f}. Meets the > $1,000,000 test: {ace_is_ai}")

    print("  - All limited partners (Liam, Jack, Ace) and the general partner (owned by them) are accredited investors.")
    print("Conclusion for A: The partnership QUALIFIES as an Accredited Investor under the look-through provision.\n")

    # --- Option B: Individual (Joint Income Test) ---
    print("B: Checking Individual (Joint Income Test)")
    joint_income_2019 = 150_000 + 170_000
    meets_2019 = joint_income_2019 > 300_000
    print(f"  - 2019 Joint Income: ${150_000:,.2f} + ${170_000:,.2f} = ${joint_income_2019:,.2f}. Meets the > $300,000 test: {meets_2019}")

    joint_income_2020 = 175_000 + 175_000
    meets_2020 = joint_income_2020 > 300_000
    print(f"  - 2020 Joint Income: ${175_000:,.2f} + ${175_000:,.2f} = ${joint_income_2020:,.2f}. Meets the > $300,000 test: {meets_2020}")
    print("  - The couple also reasonably expects to exceed this income level in the current year.")
    print("Conclusion for B: The individual QUALIFIES as an Accredited Investor.\n")

    # --- Option C: Individual (Net Asset Test) ---
    print("C: Checking Individual (Net Asset Test)")
    joint_assets = 6_000_000
    joint_liabilities = 1_000_000
    joint_net_assets = joint_assets - joint_liabilities
    meets_net_asset_test = joint_net_assets >= 5_000_000
    print(f"  - Joint Net Assets: ${joint_assets:,.2f} - ${joint_liabilities:,.2f} = ${joint_net_assets:,.2f}. Meets the >= $5,000,000 test: {meets_net_asset_test}")
    print("Conclusion for C: The individual QUALIFIES as an Accredited Investor.\n")

    # --- Option D: Corporation (Net Asset Test) ---
    print("D: Checking Corporation")
    jose_nfa = 100_000_000
    corp_assets = jose_nfa * 0.10
    corp_meets_asset_test = corp_assets >= 5_000_000
    print(f"  - Corporation's Net Assets from Jose's transfer: 10% of ${jose_nfa:,.2f} = ${corp_assets:,.2f}.")
    print(f"  - The corporation's net assets meet the >= $5,000,000 test: {corp_meets_asset_test}")
    print("  - Since the corporation qualifies on its own financial standing, it is an accredited investor.")
    print("Conclusion for D: The corporation QUALIFIES as an Accredited Investor.\n")

    # --- Option E: Corporation (Analysis of Both Tests) ---
    print("E: Checking Corporation")
    corp_e_assets = 5_500_000
    corp_e_liabilities = 300_000
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    corp_e_meets_asset_test = corp_e_net_assets >= 5_000_000
    print(f"  - Corporation's Net Assets: ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}.")
    print(f"  - On a literal basis, the corporation meets the >= $5,000,000 net asset test: {corp_e_meets_asset_test}")
    print("\n  However, let's analyze its owners (the Look-through Provision):")
    alex_is_ai = 900_000 > 1_000_000 or 3_000_000 >= 5_000_000
    print(f"  - Shareholder Alex has <$1M net financial assets and <$5M net assets. Is Alex an AI? {alex_is_ai}")
    bob_is_ai = 75_000 > 1_000_000 or 0 >= 5_000_000
    print(f"  - Shareholder Bob has <$1M net financial assets and net assets <= $0. Is Bob an AI? {bob_is_ai}")
    print("\n  - All owners of the corporation are NOT accredited investors.")
    print("  - The intent of securities law is to protect investors like Alex and Bob. Allowing a corporation, owned entirely by non-accredited individuals, to qualify simply because it holds sufficient assets would defeat the purpose of the regulation.")
    print("Conclusion for E: The corporation would likely NOT be classified as an Accredited Investor.\n")

    print("--- Final Answer ---")
    print("Based on the analysis, the entity in option E would not be classified as an Accredited Investor because it is wholly owned by non-accredited individuals and appears to be an investment vehicle, which goes against the spirit of investor protection laws.")

solve_accredited_investor_case()
<<<E>>>