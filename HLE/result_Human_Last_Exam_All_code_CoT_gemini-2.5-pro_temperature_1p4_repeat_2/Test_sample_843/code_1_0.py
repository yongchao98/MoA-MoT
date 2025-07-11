def solve_accredited_investor():
    """
    Analyzes each option to determine which would not be classified as an Accredited Investor
    under Ontario securities laws as of January 2021.
    """
    # Define the key thresholds from NI 45-106
    FINANCIAL_ASSETS_THRESHOLD = 1_000_000
    NET_ASSETS_THRESHOLD = 5_000_000
    JOINT_INCOME_THRESHOLD = 300_000
    ENTITY_NET_ASSETS_THRESHOLD = 5_000_000

    print("Analyzing each option based on the Accredited Investor criteria...\n")

    # --- Analysis of Option A: Limited Partnership ---
    print("--- Analysis of Option A (Limited Partnership) ---")
    print("An entity like a limited partnership can qualify if all its beneficial owners are accredited investors.")
    
    # Check Liam
    liam_fin_assets = 5_400_000.00
    is_liam_ai = liam_fin_assets > FINANCIAL_ASSETS_THRESHOLD
    print(f"Liam's net financial assets of ${liam_fin_assets:,.2f} are greater than the ${FINANCIAL_ASSETS_THRESHOLD:,.2f} threshold. Liam is an AI: {is_liam_ai}")

    # Check Jack
    jack_assets = 18_000_000.00
    jack_liabilities = 5_000_000.00
    jack_net_assets = jack_assets - jack_liabilities
    is_jack_ai = jack_net_assets >= NET_ASSETS_THRESHOLD
    print(f"Jack's net assets (${jack_assets:,.2f} - ${jack_liabilities:,.2f} = ${jack_net_assets:,.2f}) meet the ${NET_ASSETS_THRESHOLD:,.2f} threshold. Jack is an AI: {is_jack_ai}")

    # Check Ace
    ace_fin_assets = 25_000_000.00
    is_ace_ai = ace_fin_assets > FINANCIAL_ASSETS_THRESHOLD
    print(f"Ace's net financial assets of ${ace_fin_assets:,.2f} are greater than the ${FINANCIAL_ASSETS_THRESHOLD:,.2f} threshold. Ace is an AI: {is_ace_ai}")

    # Check the General Partner (GP) corporation
    gp_gift_per_partner = 2_000_000.00 # Correcting typo in question from 2,000,0000
    gp_net_assets = gp_gift_per_partner * 3
    is_gp_ai = gp_net_assets >= ENTITY_NET_ASSETS_THRESHOLD
    print(f"The GP's net assets (${gp_net_assets:,.2f}) meet the ${ENTITY_NET_ASSETS_THRESHOLD:,.2f} entity threshold. The GP is an AI: {is_gp_ai}")
    
    print("Result for A: Since all limited and general partners are Accredited Investors, the Limited Partnership qualifies. Option A is an Accredited Investor.\n")

    # --- Analysis of Option B: Individual (Joint Income) ---
    print("--- Analysis of Option B (Individual - Joint Income Test) ---")
    print(f"The test requires combined pre-tax income to exceed ${JOINT_INCOME_THRESHOLD:,.2f} in each of the two most recent years (2019, 2020).")
    
    individual_income_2019 = 150_000.00
    spouse_income_2019 = 170_000.00
    combined_2019 = individual_income_2019 + spouse_income_2019
    is_2019_passed = combined_2019 > JOINT_INCOME_THRESHOLD
    print(f"2019 Combined Income: ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${combined_2019:,.2f}. Test passed: {is_2019_passed}")

    individual_income_2020 = 175_000.00
    spouse_income_2020 = 175_000.00
    combined_2020 = individual_income_2020 + spouse_income_2020
    is_2020_passed = combined_2020 > JOINT_INCOME_THRESHOLD
    print(f"2020 Combined Income: ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${combined_2020:,.2f}. Test passed: {is_2020_passed}")
    
    print("Result for B: The income test is passed for both years. Option B is an Accredited Investor.\n")

    # --- Analysis of Option C: Individual (Net Assets) ---
    print("--- Analysis of Option C (Individual - Net Asset Test) ---")
    print(f"The test requires net assets (alone or with a spouse) of at least ${NET_ASSETS_THRESHOLD:,.2f}.")
    total_assets = 6_000_000.00
    total_liabilities = 1_000_000.00
    net_assets = total_assets - total_liabilities
    is_net_asset_passed = net_assets >= NET_ASSETS_THRESHOLD
    print(f"Joint Net Assets: ${total_assets:,.2f} - ${total_liabilities:,.2f} = ${net_assets:,.2f}. Test passed: {is_net_asset_passed}")
    print("Result for C: The net asset test is met. Option C is an Accredited Investor.\n")

    # --- Analysis of Option D: Corporation (Entity Net Assets) ---
    print("--- Analysis of Option D (Corporation) ---")
    print(f"A corporation can qualify on its own if its net assets are at least ${ENTITY_NET_ASSETS_THRESHOLD:,.2f}.")
    jose_nfa = 100_000_000.00
    corp_transfer_pct = 0.10
    corp_net_assets = jose_nfa * corp_transfer_pct
    is_corp_ai = corp_net_assets >= ENTITY_NET_ASSETS_THRESHOLD
    print(f"Corporation's net assets: {corp_transfer_pct:.0%} of ${jose_nfa:,.2f} = ${corp_net_assets:,.2f}. Test passed: {is_corp_ai}")
    print("Result for D: The corporation meets the entity net asset test. Option D is an Accredited Investor.\n")

    # --- Analysis of Option E: Corporation (Investment Fund Exception) ---
    print("--- Analysis of Option E (Corporation) ---")
    corp_gic = 5_500_000.00
    corp_liabilities = 300_000.00
    corp_net_assets = corp_gic - corp_liabilities
    passes_net_asset_test = corp_net_assets >= ENTITY_NET_ASSETS_THRESHOLD
    print(f"Corporation's minimum net assets: ${corp_gic:,.2f} - ${corp_liabilities:,.2f} = ${corp_net_assets:,.2f}")
    print(f"On its face, the corporation's net assets of ${corp_net_assets:,.2f} meet the ${ENTITY_NET_ASSETS_THRESHOLD:,.2f} threshold: {passes_net_asset_test}")
    print("However, the rules state that the $5,000,000 net asset test for entities does not apply if the entity is an 'investment fund'.")
    print("This corporation is owned entirely by individuals who are NOT accredited investors and its assets are composed solely of passive investments (GIC, investment properties).")
    print("This structure is characteristic of an 'investment fund', created to pool assets. As an investment fund, it cannot use the $5M net asset test and does not meet other criteria (e.g., all owners being AIs).")
    print("Result for E: Because it is effectively an investment fund comprised of non-accredited investors, it would NOT be classified as an Accredited Investor.\n")

solve_accredited_investor()
<<<E>>>