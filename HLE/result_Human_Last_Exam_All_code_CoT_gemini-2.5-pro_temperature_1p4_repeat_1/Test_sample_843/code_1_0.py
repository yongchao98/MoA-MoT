import locale

# Use locale to format numbers with commas for readability
locale.setlocale(locale.LC_ALL, '')

def format_currency(value):
    return locale.format_string("%d", value, grouping=True)

def analyze_accredited_investors():
    """
    Analyzes each option based on Ontario's Accredited Investor criteria as of Jan 2021.
    """
    # Define AI thresholds
    INDIVIDUAL_FINANCIAL_ASSETS = 1000000
    INDIVIDUAL_NET_ASSETS = 5000000
    INDIVIDUAL_INCOME = 200000
    JOINT_INCOME = 300000
    ENTITY_NET_ASSETS = 5000000

    print("--- Analysis of Accredited Investor Status ---\n")

    # --- Option A: Limited Partnership ---
    print("--- Option A: Limited Partnership ---")
    liam_nfa = 5400000
    jack_assets = 18000000
    jack_liabilities = 5000000
    ace_nfa = 25000000

    print(f"Checking limited partner Liam: Net financial assets of ${format_currency(liam_nfa)}.")
    print(f"Equation: ${format_currency(liam_nfa)} > ${format_currency(INDIVIDUAL_FINANCIAL_ASSETS)} -> Status: Accredited Investor.\n")

    jack_net_assets = jack_assets - jack_liabilities
    print(f"Checking limited partner Jack: Net assets based on total assets of ${format_currency(jack_assets)} and liabilities of ${format_currency(jack_liabilities)}.")
    print(f"Equation: ${format_currency(jack_assets)} - ${format_currency(jack_liabilities)} = ${format_currency(jack_net_assets)}.")
    print(f"Result: ${format_currency(jack_net_assets)} > ${format_currency(INDIVIDUAL_NET_ASSETS)} -> Status: Accredited Investor.\n")

    print(f"Checking limited partner Ace: Net financial assets of ${format_currency(ace_nfa)}.")
    print(f"Equation: ${format_currency(ace_nfa)} > ${format_currency(INDIVIDUAL_FINANCIAL_ASSETS)} -> Status: Accredited Investor.\n")
    
    print("Checking General Partner: All owners of the GP (Liam, Jack, Ace) are Accredited Investors.")
    print("Therefore, the General Partner is an Accredited Investor via the look-through provision.\n")

    print("Conclusion for A: Since all limited and general partners are Accredited Investors, the Limited Partnership qualifies as an Accredited Investor.")
    print("-----------------------------------------\n")


    # --- Option B: Individual with Spouse (Income Test) ---
    print("--- Option B: Individual with Spouse (Income Test) ---")
    spouse_income_2019 = 170000
    spouse_income_2020 = 175000
    individual_income_2019 = 150000
    individual_income_2020 = 175000
    
    joint_income_2019 = individual_income_2019 + spouse_income_2019
    joint_income_2020 = individual_income_2020 + spouse_income_2020

    print("Checking joint income test for the last two years.")
    print(f"2019 Joint Income Equation: ${format_currency(individual_income_2019)} (individual) + ${format_currency(spouse_income_2019)} (spouse) = ${format_currency(joint_income_2019)}.")
    print(f"Result 2019: ${format_currency(joint_income_2019)} > ${format_currency(JOINT_INCOME)}.\n")

    print(f"2020 Joint Income Equation: ${format_currency(individual_income_2020)} (individual) + ${format_currency(spouse_income_2020)} (spouse) = ${format_currency(joint_income_2020)}.")
    print(f"Result 2020: ${format_currency(joint_income_2020)} > ${format_currency(JOINT_INCOME)}.\n")

    print("Conclusion for B: The individual's joint income with their spouse exceeded $300,000 in the last two years. The individual qualifies as an Accredited Investor.")
    print("-----------------------------------------\n")


    # --- Option C: Individual with Spouse (Net Assets Test) ---
    print("--- Option C: Individual with Spouse (Net Assets Test) ---")
    total_assets = 6000000
    total_liabilities = 1000000
    net_assets = total_assets - total_liabilities

    print(f"Checking joint net assets test based on total assets of ${format_currency(total_assets)} and total liabilities of ${format_currency(total_liabilities)}.")
    print(f"Equation: ${format_currency(total_assets)} - ${format_currency(total_liabilities)} = ${format_currency(net_assets)}.")
    print(f"Result: ${format_currency(net_assets)} >= ${format_currency(INDIVIDUAL_NET_ASSETS)}.\n")

    print("Conclusion for C: The individual's joint net assets with their spouse are at least $5,000,000. The individual qualifies as an Accredited Investor.")
    print("-----------------------------------------\n")


    # --- Option D: Corporation with Wealthy Owner ---
    print("--- Option D: Corporation with Wealthy Owner ---")
    jose_nfa = 100000000
    transfer_pct = 0.10
    corp_net_assets = jose_nfa * transfer_pct

    print("Checking the corporation's net assets test.")
    print(f"Jose transferred 10% of his net financial assets of ${format_currency(jose_nfa)} to the corporation.")
    print(f"Equation: ${format_currency(jose_nfa)} * {transfer_pct} = ${format_currency(corp_net_assets)}.")
    print(f"Result: The corporation's net assets are ${format_currency(corp_net_assets)}, which is greater than the ${format_currency(ENTITY_NET_ASSETS)} threshold.\n")
    
    print("Conclusion for D: The corporation qualifies as an Accredited Investor based on its own net assets, even though not all its shareholders are accredited investors.")
    print("-----------------------------------------\n")


    # --- Option E: Corporation with Ambiguous Financials ---
    print("--- Option E: Corporation with Ambiguous Financials ---")
    print("First, checking the look-through provision by examining the owners:")
    print("Owner Alex: Net financial assets of $900,000 are less than $1,000,000. Net assets of $3,000,000 are less than $5,000,000. -> Not an Accredited Investor.")
    print("Owner Bob: Insufficient income and assets. -> Not an Accredited Investor.")
    print("Result: Since not all owners are Accredited Investors, the corporation does not qualify via the look-through provision.\n")

    print("Second, checking the corporation's own net assets test (net assets > $5,000,000).")
    gic_asset = 5500000
    fin_liabilities = 300000
    print("The corporation's 'net assets' are its TOTAL assets minus its TOTAL liabilities.")
    print("Known Assets: Cashable GIC of ${:n}".format(gic_asset))
    print("Known Liabilities: Financial liabilities of ${:n}".format(fin_liabilities))
    print("The corporation also owns 'investment properties' of UNDETERMINED value, and may have further liabilities (e.g., mortgages) associated with them that are not specified.\n")

    print("The value that can be calculated is net *financial* assets, not total net assets.")
    net_financial_assets = gic_asset - fin_liabilities
    print(f"Net Financial Assets Equation: ${format_currency(gic_asset)} - ${format_currency(fin_liabilities)} = ${format_currency(net_financial_assets)}.")
    
    print("\nCrucial Point: The accredited investor test for a corporation is based on TOTAL net assets, not just financial assets. Because the value of the investment properties (a non-financial asset) and their associated liabilities are unknown, we cannot confirm that the corporation's TOTAL net assets are at least $5,000,000.")
    print("\nConclusion for E: Due to the incomplete financial information, the corporation cannot prove it meets the $5,000,000 net asset threshold. Therefore, it would NOT be classified as an Accredited Investor.")
    print("-----------------------------------------\n")

    print("Final Answer: Based on the analysis, the entity in Option E would not be classified as an Accredited Investor.")

if __name__ == "__main__":
    analyze_accredited_investors()
    print("\n<<<E>>>")