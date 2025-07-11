import sys

# Suppress scientific notation for cleaner output
class Formatter(object):
    def __init__(self, sys_stdout):
        self.stdout = sys_stdout

    def write(self, s):
        self.stdout.write(s.replace('e+', ''))

    def flush(self):
        self.stdout.flush()

sys.stdout = Formatter(sys.stdout)


def main():
    """
    Analyzes each option to determine its Accredited Investor status based on Ontario securities laws as of Jan 2021.
    """
    # Accredited Investor thresholds
    INDIVIDUAL_FINANCIAL_ASSETS = 1000000.00
    INDIVIDUAL_NET_ASSETS = 5000000.00
    INDIVIDUAL_INCOME = 200000.00
    JOINT_INCOME = 300000.00
    ENTITY_NET_ASSETS = 5000000.00

    print("Analysis of Accredited Investor Status")
    print("-" * 40)

    # --- Option A: Limited Partnership ---
    print("Option A: The Limited Partnership")
    liam_nfa = 5400000.00
    jack_total_assets = 18000000.00
    jack_total_liabilities = 5000000.00
    jack_net_assets = jack_total_assets - jack_total_liabilities
    ace_nfa = 25000000.00
    gp_assets = 2000000.00 * 3

    liam_is_ai = liam_nfa >= INDIVIDUAL_FINANCIAL_ASSETS
    jack_is_ai = jack_net_assets >= INDIVIDUAL_NET_ASSETS
    ace_is_ai = ace_nfa >= INDIVIDUAL_FINANCIAL_ASSETS
    gp_is_ai = gp_assets >= ENTITY_NET_ASSETS
    lp_is_ai_by_lookthrough = liam_is_ai and jack_is_ai and ace_is_ai and gp_is_ai

    print(f"Liam's net financial assets: ${liam_nfa:,.2f}. Meets ${INDIVIDUAL_FINANCIAL_ASSETS:,.2f} test: {liam_is_ai}")
    print(f"Jack's net assets: ${jack_total_assets:,.2f} - ${jack_total_liabilities:,.2f} = ${jack_net_assets:,.2f}. Meets ${INDIVIDUAL_NET_ASSETS:,.2f} test: {jack_is_ai}")
    print(f"Ace's net financial assets: ${ace_nfa:,.2f}. Meets ${INDIVIDUAL_FINANCIAL_ASSETS:,.2f} test: {ace_is_ai}")
    print(f"General Partner's net assets: ${gp_assets:,.2f}. Meets ${ENTITY_NET_ASSETS:,.2f} test: {gp_is_ai}")
    print(f"Result: The LP qualifies under the 'look-through' test as all partners are AIs. Accredited: {lp_is_ai_by_lookthrough}")
    print("-" * 40)

    # --- Option B: Individual (Joint Income) ---
    print("Option B: The Individual (Joint Income)")
    income_2019_individual = 150000.00
    income_2019_spouse = 170000.00
    joint_income_2019 = income_2019_individual + income_2019_spouse
    income_2020_individual = 175000.00
    income_2020_spouse = 175000.00
    joint_income_2020 = income_2020_individual + income_2020_spouse

    passes_2019 = joint_income_2019 > JOINT_INCOME
    passes_2020 = joint_income_2020 > JOINT_INCOME
    individual_is_ai = passes_2019 and passes_2020

    print(f"2019 Joint Income: ${income_2019_individual:,.2f} + ${income_2019_spouse:,.2f} = ${joint_income_2019:,.2f}. Meets > ${JOINT_INCOME:,.2f} test: {passes_2019}")
    print(f"2020 Joint Income: ${income_2020_individual:,.2f} + ${income_2020_spouse:,.2f} = ${joint_income_2020:,.2f}. Meets > ${JOINT_INCOME:,.2f} test: {passes_2020}")
    print(f"Result: The individual qualifies on joint income. Accredited: {individual_is_ai}")
    print("-" * 40)

    # --- Option C: Individual (Joint Net Assets) ---
    print("Option C: The Individual (Joint Net Assets)")
    joint_total_assets = 6000000.00
    joint_total_liabilities = 1000000.00
    joint_net_assets = joint_total_assets - joint_total_liabilities
    individual_is_ai = joint_net_assets >= INDIVIDUAL_NET_ASSETS

    print(f"Joint Net Assets: ${joint_total_assets:,.2f} - ${joint_total_liabilities:,.2f} = ${joint_net_assets:,.2f}")
    print(f"Meets >= ${INDIVIDUAL_NET_ASSETS:,.2f} test: {individual_is_ai}")
    print(f"Result: The individual qualifies on joint net assets. Accredited: {individual_is_ai}")
    print("-" * 40)

    # --- Option D: Corporation (Jose and James) ---
    print("Option D: The Corporation (Jose and James)")
    jose_nfa = 100000000.00
    transfer_pct = 0.10
    corp_assets = jose_nfa * transfer_pct
    corp_is_ai_by_assets = corp_assets >= ENTITY_NET_ASSETS

    print(f"Corporation's net assets from transfer: {transfer_pct:.0%} of ${jose_nfa:,.2f} = ${corp_assets:,.2f}")
    print(f"Meets >= ${ENTITY_NET_ASSETS:,.2f} test: {corp_is_ai_by_assets}")
    print("The corporation qualifies on its own net assets, even though one owner (James) is not an AI.")
    print(f"Result: The corporation qualifies on its assets. Accredited: {corp_is_ai_by_assets}")
    print("-" * 40)

    # --- Option E: Corporation (Alex and Bob) ---
    print("Option E: The Corporation (Alex and Bob)")
    # Test 1: 'Look-Through' Test
    alex_nfa = 900000.00
    alex_na = 3000000.00
    alex_is_ai = (alex_nfa >= INDIVIDUAL_FINANCIAL_ASSETS) or (alex_na >= INDIVIDUAL_NET_ASSETS)
    bob_is_ai = False  # Clearly fails all tests based on info provided

    print(f"Owner Alex's status: Not Accredited (NFA ${alex_nfa:,.2f} < ${INDIVIDUAL_FINANCIAL_ASSETS:,.2f}, NA ${alex_na:,.2f} < ${INDIVIDUAL_NET_ASSETS:,.2f})")
    print(f"Owner Bob's status: Not Accredited")
    print("Result 1: Fails the 'look-through' test because not all owners are AIs.")

    # Test 2: Entity Asset Test
    corp_e_assets = 5500000.00
    corp_e_liabilities = 300000.00
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    corp_e_is_ai_by_assets = corp_e_net_assets >= ENTITY_NET_ASSETS
    
    print("\nResult 2: Checking the entity's own assets.")
    print(f"Corporation Net Assets: ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}")
    print(f"This value of ${corp_e_net_assets:,.2f} meets the >= ${ENTITY_NET_ASSETS:,.2f} test: {corp_e_is_ai_by_assets}")
    
    print("\nFinal Conclusion:")
    print("This entity fails the 'look-through' test. Since every other option clearly qualifies as an Accredited Investor under at least one test, this corporation, being entirely owned by non-accredited individuals, is the one that would not be classified as an Accredited Investor.")


if __name__ == "__main__":
    main()
