import locale

# Use locale to format numbers with commas for readability
try:
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
except locale.Error:
    try:
        locale.setlocale(locale.LC_ALL, 'English_United States.1252')
    except locale.Error:
        pass # Locale settings not available, use default

def format_currency(value):
    """Formats a number as a currency string."""
    return locale.currency(value, symbol=True, grouping=True)

def print_result(description, num1, operation, num2, comparison, threshold, passed):
    """Helper function to print calculation results."""
    status = "Pass" if passed else "Fail"
    print(f"- {description}: {format_currency(num1)} {operation} {format_currency(num2)} = {format_currency(num1 + num2 if operation == '+' else num1 - num2)} {comparison} {format_currency(threshold)} [{status}]")

# --- Accredited Investor Thresholds ---
INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD = 1_000_000
INDIVIDUAL_NET_ASSET_THRESHOLD = 5_000_000
JOINT_INCOME_THRESHOLD = 300_000
ENTITY_NET_ASSET_THRESHOLD = 5_000_000

def analyze_options():
    """Analyzes each option to determine its accredited investor status."""
    
    print("Analyzing each option to determine accredited investor status...\n")

    # --- Option A: Limited Partnership ---
    print("--- Option A: Limited Partnership ---")
    print("A limited partnership can qualify if all its limited partners are accredited investors.")
    # Liam's Status
    liam_nfa = 5_400_000
    liam_is_ai = liam_nfa > INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD
    print(f"- Liam's Net Financial Assets: {format_currency(liam_nfa)} > {format_currency(INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD)} [{ 'Pass' if liam_is_ai else 'Fail'}]")
    
    # Jack's Status
    jack_assets = 18_000_000
    jack_liabilities = 5_000_000
    jack_net_assets = jack_assets - jack_liabilities
    jack_is_ai = jack_net_assets >= INDIVIDUAL_NET_ASSET_THRESHOLD
    print_result("Jack's Net Assets", jack_assets, '-', jack_liabilities, '>=', INDIVIDUAL_NET_ASSET_THRESHOLD, jack_is_ai)

    # Ace's Status
    ace_nfa = 25_000_000
    ace_is_ai = ace_nfa > INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD
    print(f"- Ace's Net Financial Assets: {format_currency(ace_nfa)} > {format_currency(INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD)} [{ 'Pass' if ace_is_ai else 'Fail'}]")
    
    is_a_ai = liam_is_ai and jack_is_ai and ace_is_ai
    print(f"Conclusion: {'All limited partners are accredited investors, so the partnership QUALIFIES.' if is_a_ai else 'NOT all limited partners are AIs, so the partnership does NOT QUALIFY.'}\n")

    # --- Option B: Individual (Joint Income) ---
    print("--- Option B: Individual (Joint Income) ---")
    print("An individual can qualify if their joint income with a spouse exceeded $300,000 in each of the last two years.")
    income_2019_ind = 150_000
    income_2019_spouse = 170_000
    joint_income_2019 = income_2019_ind + income_2019_spouse
    passed_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
    print_result("Joint Income 2019", income_2019_ind, '+', income_2019_spouse, '>', JOINT_INCOME_THRESHOLD, passed_2019)
    
    income_2020_ind = 175_000
    income_2020_spouse = 175_000
    joint_income_2020 = income_2020_ind + income_2020_spouse
    passed_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD
    print_result("Joint Income 2020", income_2020_ind, '+', income_2020_spouse, '>', JOINT_INCOME_THRESHOLD, passed_2020)
    
    is_b_ai = passed_2019 and passed_2020
    print(f"Conclusion: {'The joint income test is met for both years, so the individual QUALIFIES.' if is_b_ai else 'The joint income test is not met.'}\n")

    # --- Option C: Individual (Joint Net Assets) ---
    print("--- Option C: Individual (Joint Net Assets) ---")
    print("An individual can qualify if their joint net assets with a spouse are at least $5,000,000.")
    assets = 6_000_000
    liabilities = 1_000_000
    net_assets = assets - liabilities
    is_c_ai = net_assets >= INDIVIDUAL_NET_ASSET_THRESHOLD
    print_result("Joint Net Assets", assets, '-', liabilities, '>=', INDIVIDUAL_NET_ASSET_THRESHOLD, is_c_ai)
    print(f"Conclusion: {'The joint net asset test is met, so the individual QUALIFIES.' if is_c_ai else 'The joint net asset test is NOT met.'}\n")
    
    # --- Option D: Corporation (Net Assets/Ownership) ---
    print("--- Option D: Corporation ---")
    print("A corporation can qualify on its own net assets (>= $5M) or if a majority of its shares are owned by accredited investors.")
    jose_nfa = 100_000_000
    corp_assets_from_jose = jose_nfa * 0.10
    corp_is_ai_on_assets = corp_assets_from_jose >= ENTITY_NET_ASSET_THRESHOLD
    print(f"- Corporation's Net Assets from Jose's transfer ({format_currency(jose_nfa)} * 10%): {format_currency(corp_assets_from_jose)} >= {format_currency(ENTITY_NET_ASSET_THRESHOLD)} [{ 'Pass' if corp_is_ai_on_assets else 'Fail'}]")

    jose_is_ai = jose_nfa > INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD
    print(f"- Majority shareholder Jose's Net Financial Assets: {format_currency(jose_nfa)} > {format_currency(INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD)} [{ 'Pass' if jose_is_ai else 'Fail'}]")
    print(f"Conclusion: {'The corporation qualifies on its own net assets AND because its majority shareholder (99%) is an AI. It QUALIFIES.' if corp_is_ai_on_assets and jose_is_ai else 'The corporation does NOT qualify.'}\n")

    # --- Option E: Corporation (Net Assets & Solely Created Rule) ---
    print("--- Option E: Corporation ---")
    print("This corporation must be assessed on its net assets and the 'solely created' anti-avoidance rule.")
    
    alex_nfa = 900_000
    alex_na = 3_000_000
    alex_is_ai = alex_nfa > INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD or alex_na >= INDIVIDUAL_NET_ASSET_THRESHOLD
    print(f"- Owner Alex's Net Financial Assets: {format_currency(alex_nfa)} is less than {format_currency(INDIVIDUAL_NET_FINANCIAL_ASSET_THRESHOLD)}. Alex is not an AI.")
    
    print("- Owner Bob's assets and income are below all thresholds. Bob is not an AI.")
    print("Conclusion on ownership: Since neither owner is an AI, the corporation cannot qualify based on its ownership.")
    
    corp_assets = 5_500_000
    corp_liabilities = 300_000
    corp_net_assets = corp_assets - corp_liabilities
    passes_asset_test = corp_net_assets >= ENTITY_NET_ASSET_THRESHOLD
    print_result("\nCorporation's Net Assets (GIC minus liabilities)", corp_assets, '-', corp_liabilities, '>=', ENTITY_NET_ASSET_THRESHOLD, passes_asset_test)
    
    print("\nHowever, the 'solely created' rule must be considered.")
    print("This rule can disqualify a corporation that qualifies ONLY on its net assets if it was created or is used solely to purchase securities as an AI.")
    print("Given that the owners are not AIs and the corporation appears to be a passive investment vehicle, it is likely disqualified under this rule.")
    
    print(f"Conclusion: The corporation is structured to meet the net asset test, but because its owners are not AIs and its purpose appears to be passive investment, it would NOT be classified as an Accredited Investor.\n")
    
    final_answer = "E"
    print(f"The option that would NOT be classified as an Accredited Investor is E.")
    return final_answer
    
if __name__ == '__main__':
    final_answer = analyze_options()
    # The final answer format is specific
    # print(f"<<<{final_answer}>>>")
    
