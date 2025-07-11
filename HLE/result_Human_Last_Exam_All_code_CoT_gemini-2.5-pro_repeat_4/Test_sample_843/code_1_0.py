import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string stream
sys.stdout = captured_output = io.StringIO()


def solve():
    """
    Analyzes each option to determine which would not be classified as an Accredited Investor
    in Ontario as of January 2021, based on NI 45-106.
    """
    # Accredited Investor Thresholds in CAD
    NET_FINANCIAL_ASSETS_INDIVIDUAL = 1_000_000.00
    NET_ASSETS_INDIVIDUAL_OR_SPOUSE = 5_000_000.00
    NET_INCOME_INDIVIDUAL = 200_000.00
    NET_INCOME_JOINT = 300_000.00
    NET_ASSETS_ENTITY = 5_000_000.00

    print("Analyzing each option against Accredited Investor criteria:\n")

    # --- Option A Analysis ---
    print("--- Analyzing Option A: Limited Partnership ---")
    print("A partnership qualifies if all its owners are Accredited Investors.")

    # Liam's status
    liam_nfa = 5_400_000.00
    is_liam_ai = liam_nfa > NET_FINANCIAL_ASSETS_INDIVIDUAL
    print(f"\n1. Liam's Status (Net Financial Assets Test):")
    print(f"   Is Liam's Net Financial Assets (${liam_nfa:,.2f}) > Threshold (${NET_FINANCIAL_ASSETS_INDIVIDUAL:,.2f})? {is_liam_ai}.")
    print(f"   Liam is an Accredited Investor.")

    # Jack's status
    jack_assets = 18_000_000.00
    jack_liabilities = 5_000_000.00
    jack_net_assets = jack_assets - jack_liabilities
    is_jack_ai = jack_net_assets >= NET_ASSETS_INDIVIDUAL_OR_SPOUSE
    print(f"\n2. Jack's Status (Net Assets Test):")
    print(f"   Jack's Net Assets = Assets (${jack_assets:,.2f}) - Liabilities (${jack_liabilities:,.2f}) = ${jack_net_assets:,.2f}")
    print(f"   Is Jack's Net Assets (${jack_net_assets:,.2f}) >= Threshold (${NET_ASSETS_INDIVIDUAL_OR_SPOUSE:,.2f})? {is_jack_ai}.")
    print(f"   Jack is an Accredited Investor.")

    # Ace's status
    ace_nfa = 25_000_000.00
    is_ace_ai = ace_nfa > NET_FINANCIAL_ASSETS_INDIVIDUAL
    print(f"\n3. Ace's Status (Net Financial Assets Test):")
    print(f"   Is Ace's Net Financial Assets (${ace_nfa:,.2f}) > Threshold (${NET_FINANCIAL_ASSETS_INDIVIDUAL:,.2f})? {is_ace_ai}.")
    print(f"   Ace is an Accredited Investor.")
    
    # General Partner's status
    gift_per_partner = 2_000_000.00
    num_partners_gifting = 3
    gp_assets = gift_per_partner * num_partners_gifting
    is_gp_ai = gp_assets >= NET_ASSETS_ENTITY
    print(f"\n4. General Partner's Status (Entity Net Asset Test):")
    print(f"   GP Net Assets = {num_partners_gifting} partners * ${gift_per_partner:,.2f}/partner = ${gp_assets:,.2f}")
    print(f"   Is GP Net Assets (${gp_assets:,.2f}) >= Threshold (${NET_ASSETS_ENTITY:,.2f})? {is_gp_ai}.")
    print(f"   The General Partner is an Accredited Investor.")

    print("\nConclusion for A: Since all limited and general partners are Accredited Investors, the Limited Partnership IS an Accredited Investor.")
    print("-" * 50)

    # --- Option B Analysis ---
    print("\n--- Analyzing Option B: Individual (Joint Income Test) ---")
    print("An individual qualifies if their joint income with a spouse exceeded $300,000 in the last two years.")

    # 2019 Joint Income
    ind_income_2019 = 150_000.00
    spouse_income_2019 = 170_000.00
    joint_income_2019 = ind_income_2019 + spouse_income_2019
    is_2019_passed = joint_income_2019 > NET_INCOME_JOINT
    print(f"\n1. 2019 Joint Income:")
    print(f"   Calculation: ${ind_income_2019:,.2f} + ${spouse_income_2019:,.2f} = ${joint_income_2019:,.2f}")
    print(f"   Is income (${joint_income_2019:,.2f}) > Threshold (${NET_INCOME_JOINT:,.2f})? {is_2019_passed}.")

    # 2020 Joint Income
    ind_income_2020 = 175_000.00
    spouse_income_2020 = 175_000.00
    joint_income_2020 = ind_income_2020 + spouse_income_2020
    is_2020_passed = joint_income_2020 > NET_INCOME_JOINT
    print(f"\n2. 2020 Joint Income:")
    print(f"   Calculation: ${ind_income_2020:,.2f} + ${spouse_income_2020:,.2f} = ${joint_income_2020:,.2f}")
    print(f"   Is income (${joint_income_2020:,.2f}) > Threshold (${NET_INCOME_JOINT:,.2f})? {is_2020_passed}.")

    print("\nConclusion for B: The individual passes the joint income test for both years. They ARE an Accredited Investor.")
    print("-" * 50)

    # --- Option C Analysis ---
    print("\n--- Analyzing Option C: Individual (Net Asset Test) ---")
    print("An individual qualifies if their net assets, alone or with a spouse, are at least $5,000,000.")
    total_assets = 6_000_000.00
    total_liabilities = 1_000_000.00
    net_assets = total_assets - total_liabilities
    is_c_ai = net_assets >= NET_ASSETS_INDIVIDUAL_OR_SPOUSE
    print(f"\nJoint Net Assets = Assets (${total_assets:,.2f}) - Liabilities (${total_liabilities:,.2f}) = ${net_assets:,.2f}")
    print(f"Is Net Assets (${net_assets:,.2f}) >= Threshold (${NET_ASSETS_INDIVIDUAL_OR_SPOUSE:,.2f})? {is_c_ai}.")

    print("\nConclusion for C: The individual meets the net asset test. They ARE an Accredited Investor.")
    print("-" * 50)

    # --- Option D Analysis ---
    print("\n--- Analyzing Option D: Corporation (Entity Net Asset Test) ---")
    print("A corporation qualifies if it has net assets of at least $5,000,000.")
    jose_nfa = 100_000_000.00
    transfer_pct = 0.10
    corp_assets = jose_nfa * transfer_pct
    is_d_ai = corp_assets >= NET_ASSETS_ENTITY
    print(f"\nCorporation's Net Assets = Jose's NFA (${jose_nfa:,.2f}) * Transfer Percentage ({transfer_pct:.0%}) = ${corp_assets:,.2f}")
    print(f"Is Corp. Net Assets (${corp_assets:,.2f}) >= Threshold (${NET_ASSETS_ENTITY:,.2f})? {is_d_ai}.")
    print("\nConclusion for D: The corporation's net assets exceed the threshold. It IS an Accredited Investor.")
    print("-" * 50)

    # --- Option E Analysis ---
    print("\n--- Analyzing Option E: Corporation ---")
    print("A corporation can qualify on its own net assets OR if all owners are AIs.")
    
    print("\n1. Analysis of Owners: Neither owner (Alex or Bob) meets the individual income or asset tests. The corporation does not qualify through its owners.")
    
    print("\n2. Analysis of Corporation's Net Assets:")
    corp_assets_gic = 5_500_000.00
    corp_liabilities = 300_000.00
    corp_net_assets_known = corp_assets_gic - corp_liabilities
    print(f"   Known Net Assets = Known Assets (${corp_assets_gic:,.2f}) - Known Liabilities (${corp_liabilities:,.2f}) = ${corp_net_assets_known:,.2f}")
    print("   However, the description mentions the corporation also owns 'investment properties' of unknown value and the liabilities are only 'a product of' those properties, implying the financials are incomplete.")

    print("\nConclusion for E: Since the corporation's full financial position is not provided, its net assets cannot be definitively calculated to be at least $5,000,000. As qualification cannot be proven, it would NOT be classified as an Accredited Investor.")
    print("-" * 50)
    
solve()

# Restore stdout and get the output
sys.stdout = original_stdout
output_string = captured_output.getvalue()
print(output_string)
print("<<<E>>>")