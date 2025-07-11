import sys
# Redirect stdout to a string buffer to control the final output format.
# This is a helper to ensure the final output fits the required format.
import io
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()


def analyze_investor_status():
    """
    Analyzes five different scenarios to determine which would not be classified as an Accredited Investor (AI)
    under Ontario securities laws as of January 2021.
    """

    # --- Accredited Investor Thresholds ---
    INDIVIDUAL_NET_INCOME = 200000
    JOINT_NET_INCOME = 300000
    NET_FINANCIAL_ASSETS = 1000000
    NET_ASSETS = 5000000
    ENTITY_NET_ASSETS = 5000000

    print("--- Analysis of Accredited Investor Status ---\n")

    # --- Option A: The Limited Partnership ---
    print("--- Option A: Limited Partnership Analysis ---")
    # All owners must be AIs for the LP to qualify via the look-through test.
    # Liam's status:
    liam_net_financial_assets = 5400000
    is_liam_ai = liam_net_financial_assets > NET_FINANCIAL_ASSETS
    print(f"Liam's net financial assets are ${liam_net_financial_assets:,.2f}. This is greater than the ${NET_FINANCIAL_ASSETS:,.2f} threshold. Liam is an AI.")

    # Jack's status:
    jack_total_assets = 18000000
    jack_total_liabilities = 5000000
    jack_net_assets = jack_total_assets - jack_total_liabilities
    is_jack_ai = jack_net_assets >= NET_ASSETS
    print(f"Jack's net assets are ${jack_total_assets:,.2f} - ${jack_total_liabilities:,.2f} = ${jack_net_assets:,.2f}. This is greater than the ${NET_ASSETS:,.2f} threshold. Jack is an AI.")

    # Ace's status:
    ace_net_financial_assets = 25000000
    is_ace_ai = ace_net_financial_assets > NET_FINANCIAL_ASSETS
    print(f"Ace's net financial assets are ${ace_net_financial_assets:,.2f}. This is greater than the ${NET_FINANCIAL_ASSETS:,.2f} threshold. Ace is an AI.")
    
    # General Partner's status:
    gp_assets = 2000000 + 2000000 + 2000000
    is_gp_ai = gp_assets >= ENTITY_NET_ASSETS
    print(f"The General Partner's net assets are ${gp_assets:,.2f}. This is greater than the ${ENTITY_NET_ASSETS:,.2f} threshold. The GP is an AI.")

    is_a_ai = is_liam_ai and is_jack_ai and is_ace_ai and is_gp_ai
    print(f"Conclusion: All owners of the limited partnership are Accredited Investors. Therefore, the LP is an Accredited Investor.\n")

    # --- Option B: The Individual (Joint Income) ---
    print("--- Option B: Individual Analysis (Joint Income) ---")
    individual_income_2019 = 150000
    individual_income_2020 = 175000
    spouse_income_2019 = 170000
    spouse_income_2020 = 175000

    joint_income_2019 = individual_income_2019 + spouse_income_2019
    joint_income_2020 = individual_income_2020 + spouse_income_2020
    
    print(f"Joint income in 2019: ${individual_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}")
    print(f"Joint income in 2020: ${individual_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}")
    
    is_b_ai = joint_income_2019 > JOINT_NET_INCOME and joint_income_2020 > JOINT_NET_INCOME
    print(f"Conclusion: The joint income exceeded the ${JOINT_NET_INCOME:,.2f} threshold in the two most recent years. The individual is an Accredited Investor.\n")

    # --- Option C: The Individual (Joint Net Assets) ---
    print("--- Option C: Individual Analysis (Joint Net Assets) ---")
    joint_total_assets = 6000000
    joint_total_liabilities = 1000000
    joint_net_assets = joint_total_assets - joint_total_liabilities
    
    print(f"Joint net assets: ${joint_total_assets:,.2f} - ${joint_total_liabilities:,.2f} = ${joint_net_assets:,.2f}")

    is_c_ai = joint_net_assets >= NET_ASSETS
    print(f"Conclusion: The joint net assets of ${joint_net_assets:,.2f} meet the ${NET_ASSETS:,.2f} threshold. The individual is an Accredited Investor.\n")

    # --- Option D: The Corporation (with one non-AI shareholder) ---
    print("--- Option D: Corporation Analysis ---")
    # First, check if the entity itself qualifies.
    jose_nfa = 100000000
    corp_d_assets = 0.10 * jose_nfa
    
    print(f"The corporation's net assets are 10% of ${jose_nfa:,.2f}, which is ${corp_d_assets:,.2f}.")
    is_corp_d_ai_by_assets = corp_d_assets >= ENTITY_NET_ASSETS
    
    # Check the look-through provision just for completeness, although it's not needed if the entity qualifies on its own.
    james_income_2019 = 75000
    james_spouse_income_2019 = 210000
    james_joint_income_2019 = james_income_2019 + james_spouse_income_2019
    
    james_income_2020 = 75000
    james_spouse_income_2020 = 210000
    james_joint_income_2020 = james_income_2020 + james_spouse_income_2020

    print(f"James' joint income in 2019 was ${james_income_2019:,.2f} + ${james_spouse_income_2019:,.2f} = ${james_joint_income_2019:,.2f}")
    print(f"James' joint income in 2020 was ${james_income_2020:,.2f} + ${james_spouse_income_2020:,.2f} = ${james_joint_income_2020:,.2f}")
    print(f"James does not qualify as an AI because his joint income is below the ${JOINT_NET_INCOME:,.2f} threshold.")

    print(f"Conclusion: Although not all shareholders are AIs, the corporation itself qualifies as an AI because its net assets of ${corp_d_assets:,.2f} exceed the ${ENTITY_NET_ASSETS:,.2f} threshold. The corporation is an Accredited Investor.\n")
    
    # --- Option E: The Corporation (with two non-AI shareholders) ---
    print("--- Option E: Corporation Analysis ---")
    # Step 1: Check the corporation's net assets.
    corp_e_assets = 5500000
    corp_e_liabilities = 300000
    corp_e_net_assets = corp_e_assets - corp_e_liabilities
    
    print(f"Corporation's net assets: ${corp_e_assets:,.2f} - ${corp_e_liabilities:,.2f} = ${corp_e_net_assets:,.2f}")
    is_corp_e_ai_by_assets = corp_e_net_assets >= ENTITY_NET_ASSETS
    print(f"On a purely numerical basis, the corporation's net assets of ${corp_e_net_assets:,.2f} are greater than the ${ENTITY_NET_ASSETS:,.2f} threshold.")

    # Step 2: Check the shareholder status.
    alex_nfa = 900000
    alex_net_assets = 3000000
    print(f"Alex's net financial assets (${alex_nfa:,.2f}) and net assets (${alex_net_assets:,.2f}) are below the AI thresholds. Alex is not an AI.")
    bob_income_2020 = 41000
    bob_nfa = 75000
    print(f"Bob's income (${bob_income_2020:,.2f}) and net financial assets (${bob_nfa:,.2f}) are below the AI thresholds. Bob is not an AI.")
    
    # Step 3: Apply the "artificial creation" rule.
    print("\nApplying the anti-avoidance rule (Companion Policy 45-106CP):")
    print("An entity created or used solely to pool funds from non-accredited investors to meet the $5 million net asset test does not qualify as an AI.")
    print("In this case, all shareholders (Alex and Bob) are not Accredited Investors. The corporation's primary holdings are investments (GIC, properties), and it just clears the net asset threshold.")
    print("Conclusion: This corporation fits the description of an artificial entity created to circumvent AI rules. Therefore, it would NOT be classified as an Accredited Investor.\n")
    
    print("--- Final Determination ---")
    print("Option E is the entity that would not be classified as an Accredited Investor.")
    
    return 'E'


# Run the analysis
final_answer = analyze_investor_status()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())
print(f"<<<{final_answer}>>>")