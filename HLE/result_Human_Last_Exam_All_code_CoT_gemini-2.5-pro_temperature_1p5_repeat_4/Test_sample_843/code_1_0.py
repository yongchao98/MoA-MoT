def analyze_accredited_investors():
    """
    Analyzes several scenarios to determine which would not be classified as an Accredited Investor in Ontario as of Jan 2021.
    """
    # Define Accredited Investor (AI) thresholds
    AI_FINANCIAL_ASSETS = 1000000
    AI_NET_ASSETS_INDIVIDUAL = 5000000
    AI_INCOME_JOINT = 300000
    AI_NET_ASSETS_ENTITY = 5000000

    print("Analyzing which option would not be classified as an Accredited Investor.\n")

    # --- Option E Analysis ---
    # This is the option that fails to qualify. We will show the math.
    print("--- Option E Breakdown ---")

    # Step 1: Check the status of the corporation's owners, Alex and Bob.
    print("\nStep 1: Are the owners, Alex and Bob, accredited investors?")
    alex_net_financial_assets = 900000
    alex_net_assets = 3000000
    alex_is_ai = (alex_net_financial_assets >= AI_FINANCIAL_ASSETS) or (alex_net_assets >= AI_NET_ASSETS_INDIVIDUAL)
    print(f"Alex's Net Financial Assets Test: ${alex_net_financial_assets:,.2f} < ${AI_FINANCIAL_ASSETS:,.2f} (Fail)")
    print(f"Alex's Net Asset Test: ${alex_net_assets:,.2f} < ${AI_NET_ASSETS_INDIVIDUAL:,.2f} (Fail)")
    print(f"Result: Alex is not an accredited investor.\n")
    
    bob_financial_assets = 75000
    bob_net_assets = 0
    bob_is_ai = (bob_financial_assets >= AI_FINANCIAL_ASSETS) or (bob_net_assets >= AI_NET_ASSETS_INDIVIDUAL)
    print(f"Bob's Net Financial Assets Test: ${bob_financial_assets:,.2f} < ${AI_FINANCIAL_ASSETS:,.2f} (Fail)")
    print(f"Bob's Net Asset Test: ${bob_net_assets:,.2f} < ${AI_NET_ASSETS_INDIVIDUAL:,.2f} (Fail)")
    print(f"Result: Bob is not an accredited investor.\n")

    corp_e_owners_are_ai = all([alex_is_ai, bob_is_ai])
    print("Conclusion for Step 1: Since neither Alex nor Bob are accredited investors, the corporation fails the test where 'all owners are accredited investors'.\n")

    # Step 2: Check the corporation's own net assets.
    print("Step 2: Does the corporation qualify based on its own net assets (>= $5,000,000)?")
    corp_e_assets = 5500000
    corp_e_liabilities = 300000
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    
    print(f"Corporation's Net Asset Calculation: ${corp_e_assets:,.2f} (GIC Asset) - ${corp_e_liabilities:,.2f} (Liabilities) = ${corp_e_net_assets:,.2f}")
    
    corp_e_asset_test_passes = corp_e_net_assets >= AI_NET_ASSETS_ENTITY
    print(f"Does ${corp_e_net_assets:,.2f} meet the 'at least ${AI_NET_ASSETS_ENTITY:,.2f}' test? {corp_e_asset_test_passes}")
    print("\nFinal Conclusion:")
    print("Although the corporation's net assets arithmetically exceed $5,000,000, it is the only entity entirely owned by non-accredited individuals.")
    print("This structure is the most likely to be scrutinized and rejected by regulators under anti-avoidance rules, as it appears designed to allow non-accredited individuals to invest in the exempt market.")
    print("Therefore, it is the option that would not be classified as an Accredited Investor.")


analyze_accredited_investors()