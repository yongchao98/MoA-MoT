def analyze_securities_distributions():
    """
    Analyzes several distribution scenarios based on simplified Ontario securities regulations
    as of January 2020, focusing on prospectus exemptions.
    """

    # Accredited investor financial thresholds (for an individual)
    ACCREDITED_INVESTOR_INCOME_THRESHOLD = 200000
    ACCREDITED_INVESTOR_NET_ASSET_THRESHOLD = 5000000

    print("Analyzing compliance of securities distributions with Ontario regulations.\n")

    # --- Scenario A ---
    print("--- Analysis of Scenario A ---")
    print("Distribution: Bonds by Bank of China (Canada) to retail investors.")
    # Rule: Debt securities from a Schedule I or II bank are prospectus-exempt.
    # Bank of China (Canada) is a Schedule II bank.
    is_A_compliant = True
    print("Issuer type: Schedule II Bank")
    print("Security type: Debt (Bonds)")
    print("Result: This distribution IS COMPLIANT.")
    print("Reason: The Ontario Securities Act provides a prospectus exemption for debt securities issued by a bank listed in Schedule I or II of the Bank Act (Canada).\n")

    # --- Scenario B ---
    print("--- Analysis of Scenario B ---")
    print("Distribution: Bonds by JPMorgan Chase to a large number of investors.")
    # Rule: The general exemption for banks does not typically cover Schedule III banks.
    is_B_compliant = False
    print("Issuer type: Schedule III Bank")
    print("Result: This distribution is NOT COMPLIANT without using other specific exemptions.")
    print("Reason: The broad prospectus exemption for bank-issued debt does not apply to Schedule III foreign bank branches.\n")

    # --- Scenario C ---
    print("--- Analysis of Scenario C ---")
    print("Distribution: Shares by a private issuer to an unconnected individual.")
    investor_salary = 35000
    investor_net_assets = 10000
    # Rule: Investor must meet criteria, e.g., be an accredited investor.
    is_accredited = (investor_salary >= ACCREDITED_INVESTOR_INCOME_THRESHOLD or
                     investor_net_assets >= ACCREDITED_INVESTOR_NET_ASSET_THRESHOLD)
    print("Checking 'Accredited Investor' status:")
    print(f"Investor's Salary: ${investor_salary:,}")
    print(f"Required Income Threshold: ${ACCREDITED_INVESTOR_INCOME_THRESHOLD:,}")
    print(f"Investor's Net Assets: ${investor_net_assets:,}")
    print(f"Required Net Asset Threshold: ${ACCREDITED_INVESTOR_NET_ASSET_THRESHOLD:,}")
    print("Status: Investor does not qualify as an Accredited Investor.")
    print("In addition, the investor has no connection for the 'Private Issuer' exemption to apply.")
    print("Result: This distribution is NOT COMPLIANT.\n")

    # --- Scenario D ---
    print("--- Analysis of Scenario D ---")
    print("Distribution: Shares by Fairstone Bank of Canada to retail investors.")
    # Rule: Bank exemption is for debt, not equity. Fairstone is a Schedule I bank.
    is_D_compliant = False
    print("Issuer type: Schedule I Bank")
    print("Security type: Equity (Shares)")
    print("Result: This distribution is NOT COMPLIANT.")
    print("Reason: The prospectus exemption for banks applies to their debt securities, not their equity securities (shares).\n")

    # --- Scenario E ---
    print("--- Analysis of Scenario E ---")
    print("Distribution: Shares by Caisse populaire acadienne ltée to retail investors.")
    # Rule: Ontario exemption for credit unions applies to those under the Ontario governing act.
    is_E_compliant = False
    print("Issuer type: Credit Union (likely non-Ontario based on 'acadienne' name)")
    print("Result: This distribution is NOT COMPLIANT.")
    print("Reason: The issuer appears to be from New Brunswick, so the Ontario-specific prospectus exemption for credit unions would not apply.\n")

# Execute the analysis
analyze_securities_distributions()