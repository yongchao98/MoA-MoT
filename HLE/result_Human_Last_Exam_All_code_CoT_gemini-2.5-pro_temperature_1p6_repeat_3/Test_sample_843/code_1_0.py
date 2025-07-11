import locale

# Use locale to format numbers with commas for readability
try:
    locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')
except locale.Error:
    try:
        locale.setlocale(locale.LC_ALL, 'en_US')
    except locale.Error:
        print("Locale 'en_US.UTF-8' or 'en_US' not supported. Numbers will not be formatted.")

def n(val):
    """Formats a number as a currency string."""
    if val is None:
        return "N/A"
    return f"${locale.format_string('%d', val, grouping=True)}"

# --- Accredited Investor Thresholds (Ontario, Jan 2021) ---
INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD = 1_000_000
INDIVIDUAL_NET_ASSETS_THRESHOLD = 5_000_000
INDIVIDUAL_INCOME_THRESHOLD = 200_000
JOINT_INCOME_THRESHOLD = 300_000
COMPANY_NET_ASSETS_THRESHOLD = 5_000_000

def check_option_a():
    """Evaluates the Limited Partnership."""
    print("--- Evaluating Option A: The Limited Partnership ---")
    
    # An LP is accredited if all its partners are accredited investors (Look-through provision)
    # or if the LP itself has net assets of at least $5,000,000.
    
    # Check each limited partner
    liam_net_financial_assets = 5_400_000
    liam_is_accredited = liam_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD
    print(f"Liam's check (Financial Assets): {n(liam_net_financial_assets)} > {n(INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD)}? -> {liam_is_accredited}")

    jack_assets = 18_000_000
    jack_liabilities = 5_000_000
    jack_net_assets = jack_assets - jack_liabilities
    jack_is_accredited = jack_net_assets > INDIVIDUAL_NET_ASSETS_THRESHOLD
    print(f"Jack's check (Net Assets): {n(jack_assets)} - {n(jack_liabilities)} = {n(jack_net_assets)}. Is {n(jack_net_assets)} > {n(INDIVIDUAL_NET_ASSETS_THRESHOLD)}? -> {jack_is_accredited}")

    ace_net_financial_assets = 25_000_000
    ace_is_accredited = ace_net_financial_assets > INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD
    print(f"Ace's check (Financial Assets): {n(ace_net_financial_assets)} > {n(INDIVIDUAL_FINANCIAL_ASSETS_THRESHOLD)}? -> {ace_is_accredited}")

    # Check the General Partner (a corporation)
    gp_gift = 2_000_000
    gp_num_partners = 3
    gp_assets = gp_gift * gp_num_partners
    gp_is_accredited = gp_assets >= COMPANY_NET_ASSETS_THRESHOLD
    print(f"General Partner's check (Net Assets): {gp_num_partners} * {n(gp_gift)} = {n(gp_assets)}. Is {n(gp_assets)} >= {n(COMPANY_NET_ASSETS_THRESHOLD)}? -> {gp_is_accredited}")

    is_accredited = liam_is_accredited and jack_is_accredited and ace_is_accredited and gp_is_accredited
    print(f"Conclusion for A: All partners are accredited investors. Therefore, the LP is classified as an Accredited Investor.\n")
    return "Accredited Investor"

def check_option_b():
    """Evaluates the individual based on joint income."""
    print("--- Evaluating Option B: Individual (Joint Income) ---")
    
    individual_income_2019 = 150_000
    spouse_income_2019 = 170_000
    joint_income_2019 = individual_income_2019 + spouse_income_2019
    meets_2019 = joint_income_2019 > JOINT_INCOME_THRESHOLD
    print(f"Joint Income 2019: {n(individual_income_2019)} + {n(spouse_income_2019)} = {n(joint_income_2019)}. Is {n(joint_income_2019)} > {n(JOINT_INCOME_THRESHOLD)}? -> {meets_2019}")
    
    individual_income_2020 = 175_000
    spouse_income_2020 = 175_000
    joint_income_2020 = individual_income_2020 + spouse_income_2020
    meets_2020 = joint_income_2020 > JOINT_INCOME_THRESHOLD
    print(f"Joint Income 2020: {n(individual_income_2020)} + {n(spouse_income_2020)} = {n(joint_income_2020)}. Is {n(joint_income_2020)} > {n(JOINT_INCOME_THRESHOLD)}? -> {meets_2020}")

    is_accredited = meets_2019 and meets_2020
    print(f"Conclusion for B: The joint income exceeded the threshold for the past two years. Therefore, the individual is classified as an Accredited Investor.\n")
    return "Accredited Investor"

def check_option_c():
    """Evaluates the individual based on the net asset test."""
    print("--- Evaluating Option C: Individual (Net Assets) ---")
    
    # The rule for an individual's net assets is that they must *exceed* $5,000,000.
    total_assets = 6_000_000
    total_liabilities = 1_000_000
    net_assets = total_assets - total_liabilities
    is_accredited = net_assets > INDIVIDUAL_NET_ASSETS_THRESHOLD
    print(f"Net Asset check: {n(total_assets)} - {n(total_liabilities)} = {n(net_assets)}. Is {n(net_assets)} > {n(INDIVIDUAL_NET_ASSETS_THRESHOLD)}? -> {is_accredited}")
    
    print("Conclusion for C: The individual's net assets are exactly $5,000,000, which does not *exceed* the threshold. Therefore, the individual would not be classified as an Accredited Investor based on this test.\n")
    return "Not an Accredited Investor"

def check_option_d():
    """Evaluates the corporation (Jose and James)."""
    print("--- Evaluating Option D: Corporation (Jose and James) ---")

    # A corporation is accredited if it has net assets of at least $5,000,000.
    jose_nfa = 100_000_000
    transfer_pct = 0.10
    corp_assets = jose_nfa * transfer_pct
    is_accredited = corp_assets >= COMPANY_NET_ASSETS_THRESHOLD
    print(f"Corporation's Net Asset Check: {transfer_pct:.0%} of {n(jose_nfa)} = {n(corp_assets)}. Is {n(corp_assets)} >= {n(COMPANY_NET_ASSETS_THRESHOLD)}? -> {is_accredited}")
    
    print("Conclusion for D: The corporation's net assets are at least $5,000,000. Therefore, the corporation itself is classified as an Accredited Investor (regardless of its shareholders' status).\n")
    return "Accredited Investor"

def check_option_e():
    """Evaluates the corporation (Alex and Bob)."""
    print("--- Evaluating Option E: Corporation (Alex and Bob) ---")

    # A corporation is accredited if it has net assets of at least $5,000,000.
    corp_assets = 5_500_000
    corp_liabilities = 300_000
    corp_net_assets = corp_assets - corp_liabilities
    is_accredited = corp_net_assets >= COMPANY_NET_ASSETS_THRESHOLD
    print(f"Corporation's Net Asset Check: {n(corp_assets)} - {n(corp_liabilities)} = {n(corp_net_assets)}. Is {n(corp_net_assets)} >= {n(COMPANY_NET_ASSETS_THRESHOLD)}? -> {is_accredited}")
    
    print("Conclusion for E: The corporation's net assets are at least $5,000,000. Therefore, the corporation itself is classified as an Accredited Investor (regardless of its shareholders' status).\n")
    return "Accredited Investor"

def find_non_accredited_investor():
    """Main function to run all checks and find the answer."""
    results = {
        "A": check_option_a(),
        "B": check_option_b(),
        "C": check_option_c(),
        "D": check_option_d(),
        "E": check_option_e()
    }
    
    final_answer = None
    for option, status in results.items():
        if status == "Not an Accredited Investor":
            final_answer = option
            break
            
    if final_answer:
        print(f"The option that would not be classified as an Accredited Investor is {final_answer}.")

# Run the analysis
if __name__ == "__main__":
    find_non_accredited_investor()
    print("<<<C>>>")
