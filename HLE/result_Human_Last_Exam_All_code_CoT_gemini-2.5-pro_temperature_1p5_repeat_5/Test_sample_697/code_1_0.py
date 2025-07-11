def analyze_securities_regulations():
    """
    Analyzes five scenarios against Ontario securities regulations (as of Jan 2020)
    to determine which one is compliant.
    """

    print("Analyzing each option based on prospectus exemptions in Ontario securities law...\n")

    # --- Analysis of Option C (The only one with numerical data) ---
    print("--- Analysis of Option C ---")
    investor_salary = 35000
    investor_net_assets = 10000

    # Accredited Investor Thresholds for an individual
    income_threshold = 200000
    net_financial_assets_threshold = 1000000
    net_assets_threshold = 5000000

    print("Scenario C: A private issuer distribution to an individual investor with:")
    print(f"* Salary: ${investor_salary}")
    print(f"* Net assets: ${investor_net_assets}\n")

    print("Checking if the investor qualifies as an 'Accredited Investor':")

    # Check 1: Net Income Test
    is_accredited_by_income = investor_salary > income_threshold
    print(f"1. Net Income Test: Is salary of ${investor_salary} > ${income_threshold}? Result: {is_accredited_by_income}")

    # Check 2: Net Financial Asset Test
    # The problem gives 'net assets', not 'net financial assets', but it's still far below the threshold.
    is_accredited_by_fin_assets = investor_net_assets > net_financial_assets_threshold
    print(f"2. Net Financial Assets Test: Are net assets of ${investor_net_assets} > ${net_financial_assets_threshold}? Result: {is_accredited_by_fin_assets}")

    # Check 3: Net Asset Test
    is_accredited_by_net_assets = investor_net_assets > net_assets_threshold
    print(f"3. Net Assets Test: Are net assets of ${investor_net_assets} > ${net_assets_threshold}? Result: {is_accredited_by_net_assets}\n")

    print("Conclusion for C: The investor does not meet any financial test to be an accredited investor.")
    print("The scenario also states the investor has 'no connection' to the company, ruling out other exemptions like friends and family.")
    print("Therefore, this distribution without a prospectus is NON-COMPLIANT.\n")

    # --- Analysis of Other Options (Based on Issuer Type) ---

    print("--- Analysis of Other Options ---")
    print("A. Bank of China (Canada) is a Schedule II bank under the Bank Act. The prospectus exemption for securities issued by a bank applies. This is COMPLIANT.\n")

    print("B. JPMorgan Chase operates as a Schedule III bank branch. These branches are generally restricted from raising capital from the retail public in Canada. A distribution to a 'large number of investors' would likely violate its banking license. This is NON-COMPLIANT.\n")

    print("D. Fairstone Bank of Canada is a Schedule I bank. The prospectus exemption for securities issued by a bank applies. This exemption is based on the issuer (the bank), not the investor's financial status. The fact that an investor has zero net assets is irrelevant to the legality of the distribution itself. This is COMPLIANT.\n")

    print("E. Caisse populaire acadienne ltée is a credit union. A prospectus exemption exists for credit unions distributing shares to investors who become members. This scenario is likely COMPLIANT, assuming the investors become members.\n")

    # --- Final Conclusion ---
    print("--- Final Conclusion ---")
    print("Options A, D, and E all describe activities that are exempt from prospectus requirements.")
    print("However, Option D represents the most clear-cut and textbook example of compliance. It involves a domestic (Schedule I) bank and includes details (an investor with no assets) specifically designed to test the knowledge that the bank exemption is issuer-based and does not depend on the investor's financial sophistication.")

if __name__ == "__main__":
    analyze_securities_regulations()