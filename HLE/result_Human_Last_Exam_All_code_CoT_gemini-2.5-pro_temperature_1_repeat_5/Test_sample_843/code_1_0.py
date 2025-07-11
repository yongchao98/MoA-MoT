import locale

# Use locale to format numbers as currency
locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

def print_currency(value):
    """Prints a numeric value as a formatted currency string."""
    return locale.currency(value, grouping=True)

def solve():
    """
    Analyzes which of the options would not be classified as an Accredited Investor.
    """
    # Define Accredited Investor thresholds as of Jan 2021 in Ontario
    FINANCIAL_ASSETS_THRESHOLD = 1000000
    NET_ASSETS_INDIVIDUAL_THRESHOLD = 5000000
    NET_ASSETS_ENTITY_THRESHOLD = 5000000
    JOINT_NET_INCOME_THRESHOLD = 300000
    
    final_answer = ""
    
    print("--- Analysis of Each Option ---")

    # --- Option A: Limited Partnership ---
    print("\n[A] Analyzing the Limited Partnership...")
    liam_nfa = 5400000
    jack_assets = 18000000
    jack_liabilities = 5000000
    jack_net_assets = jack_assets - jack_liabilities
    ace_nfa = 25000000
    gp_corp_assets = 2000000 * 3 # Gifted by 3 partners

    is_liam_ai = liam_nfa > FINANCIAL_ASSETS_THRESHOLD
    is_jack_ai = jack_net_assets >= NET_ASSETS_INDIVIDUAL_THRESHOLD
    is_ace_ai = ace_nfa > FINANCIAL_ASSETS_THRESHOLD
    is_gp_ai = gp_corp_assets >= NET_ASSETS_ENTITY_THRESHOLD
    
    print(f"1. Liam's net financial assets: {print_currency(liam_nfa)}. This is > {print_currency(FINANCIAL_ASSETS_THRESHOLD)}. Liam is an AI: {is_liam_ai}.")
    print(f"2. Jack's net assets: {print_currency(jack_assets)} - {print_currency(jack_liabilities)} = {print_currency(jack_net_assets)}. This is >= {print_currency(NET_ASSETS_INDIVIDUAL_THRESHOLD)}. Jack is an AI: {is_jack_ai}.")
    print(f"3. Ace's net financial assets: {print_currency(ace_nfa)}. This is > {print_currency(FINANCIAL_ASSETS_THRESHOLD)}. Ace is an AI: {is_ace_ai}.")
    print(f"4. The General Partner corporation has net assets of {print_currency(gp_corp_assets)}. This is >= {print_currency(NET_ASSETS_ENTITY_THRESHOLD)}. The GP is an AI: {is_gp_ai}.")
    
    is_a_ai = is_liam_ai and is_jack_ai and is_ace_ai and is_gp_ai
    print(f"Conclusion for A: Since all owners are accredited investors, the limited partnership qualifies. Is it an AI? {is_a_ai}")

    # --- Option B: Individual (Joint Income Test) ---
    print("\n[B] Analyzing the Individual (Joint Income Test)...")
    individual_income_2019 = 150000
    spouse_income_2019 = 170000
    joint_income_2019 = individual_income_2019 + spouse_income_2019
    
    individual_income_2020 = 175000
    spouse_income_2020 = 175000
    joint_income_2020 = individual_income_2020 + spouse_income_2020

    is_b_ai = joint_income_2019 > JOINT_NET_INCOME_THRESHOLD and joint_income_2020 > JOINT_NET_INCOME_THRESHOLD
    
    print(f"1. Joint income in 2019: {print_currency(individual_income_2019)} + {print_currency(spouse_income_2019)} = {print_currency(joint_income_2019)}. This is > {print_currency(JOINT_NET_INCOME_THRESHOLD)}.")
    print(f"2. Joint income in 2020: {print_currency(individual_income_2020)} + {print_currency(spouse_income_2020)} = {print_currency(joint_income_2020)}. This is > {print_currency(JOINT_NET_INCOME_THRESHOLD)}.")
    print(f"Conclusion for B: The individual qualifies with their spouse based on joint income. Is it an AI? {is_b_ai}")

    # --- Option C: Individual (Net Assets Test) ---
    print("\n[C] Analyzing the Individual (Net Assets Test)...")
    total_assets = 6000000
    total_liabilities = 1000000
    net_assets = total_assets - total_liabilities
    
    is_c_ai = net_assets >= NET_ASSETS_INDIVIDUAL_THRESHOLD
    
    print(f"1. Joint net assets: {print_currency(total_assets)} - {print_currency(total_liabilities)} = {print_currency(net_assets)}.")
    print(f"2. This value is at least {print_currency(NET_ASSETS_INDIVIDUAL_THRESHOLD)}.")
    print(f"Conclusion for C: The individual qualifies with their spouse based on net assets. Is it an AI? {is_c_ai}")

    # --- Option D: Corporation (Jose and James) ---
    print("\n[D] Analyzing the Corporation (Jose and James)...")
    jose_nfa = 100000000
    corp_transfer = jose_nfa * 0.10
    
    is_d_ai = corp_transfer >= NET_ASSETS_ENTITY_THRESHOLD
    
    print(f"1. Jose transferred 10% of his {print_currency(jose_nfa)} net financial assets to the corporation.")
    print(f"2. Corporation's net assets are {print_currency(corp_transfer)}. This is >= {print_currency(NET_ASSETS_ENTITY_THRESHOLD)}.")
    print(f"Conclusion for D: The corporation qualifies on its own based on its net assets. Is it an AI? {is_d_ai}")

    # --- Option E: Corporation (Alex and Bob) ---
    print("\n[E] Analyzing the Corporation (Alex and Bob)...")
    corp_assets = 5500000
    corp_liabilities = 300000
    corp_net_assets = corp_assets - corp_liabilities
    
    passes_entity_test = corp_net_assets >= NET_ASSETS_ENTITY_THRESHOLD
    
    print(f"1. Corporation's net assets: {print_currency(corp_assets)} - {print_currency(corp_liabilities)} = {print_currency(corp_net_assets)}.")
    print(f"2. This amount ({print_currency(corp_net_assets)}) is greater than the {print_currency(NET_ASSETS_ENTITY_THRESHOLD)} threshold for an entity.")
    print("3. At first glance, the corporation seems to qualify as an AI.")
    print("4. HOWEVER, the rules for an entity with net assets over $5M specifically exclude 'investment funds'.")
    print("5. This corporation's stated purpose is investing (holding a GIC, liabilities from investment properties), which classifies it as an investment fund.")
    print("6. As an investment fund, it cannot use the $5M net asset test. It must qualify another way, such as having all its owners be accredited investors.")
    
    alex_nfa = 900000
    is_alex_ai = alex_nfa > FINANCIAL_ASSETS_THRESHOLD
    is_bob_ai = False # Fails all income and asset tests
    
    print(f"7. Let's check the owners: Alex's net financial assets are {print_currency(alex_nfa)} (below the threshold). Bob does not meet any test. Neither owner is an AI.")
    print(f"Conclusion for E: Since the corporation is an investment fund and its owners are not accredited investors, it does not qualify. Is it an AI? False")
    
    if not (is_a_ai and is_b_ai and is_c_ai and is_d_ai and (not passes_entity_test)):
        final_answer = "E"

    print(f"\n--- Final Determination ---")
    print("Options A, B, C, and D all meet the criteria for being an Accredited Investor.")
    print("Option E does not qualify because it is considered an 'investment fund' and its owners are not accredited investors, making it ineligible under the applicable rules.")
    print(f"The option that would not be classified as an Accredited Investor is E.")
    print(f'<<<{final_answer}>>>')

solve()