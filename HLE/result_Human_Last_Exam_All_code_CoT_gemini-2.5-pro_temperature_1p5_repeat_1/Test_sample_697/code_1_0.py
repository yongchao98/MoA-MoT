import sys

def securities_regulation_analysis():
    """
    Analyzes several securities distribution scenarios based on Ontario
    securities regulations as of January 2020.
    """

    print("Analyzing securities distributions for compliance with Ontario regulations (as of Jan 2020)...")
    print("="*80)

    # --- Accredited Investor (AI) Thresholds for an Individual ---
    ai_income_threshold = 200000
    ai_net_financial_assets_threshold = 1000000
    ai_net_assets_threshold = 5000000

    # --- Analysis of Option C Investor ---
    investor_c_salary = 35000
    investor_c_net_assets = 10000

    is_accredited = False
    if (investor_c_salary >= ai_income_threshold or
        investor_c_net_assets >= ai_net_financial_assets_threshold or
        investor_c_net_assets >= ai_net_assets_threshold):
        is_accredited = True


    # --- Analysis of Each Option ---

    print("\n--- Analysis of Option A ---")
    print("Issuer: Bank of China (Canada)")
    print("Status: Bank of China (Canada) is a Schedule II bank under Canada's Bank Act.")
    print("Exemption: As a bank, it can rely on the 'Financial Institution Exemption' (NI 45-106, s. 2.35) to issue its own securities (like bonds) without a prospectus, even to retail investors.")
    print("Conclusion: This distribution appears compliant.")

    print("\n--- Analysis of Option B ---")
    print("Issuer: JPMorgan Chase (US Parent Entity)")
    print("Status: The US parent entity is not a 'bank' as defined by Canada's Bank Act for the purposes of this exemption.")
    print("Exemption: Cannot rely on the 'Financial Institution Exemption' for a broad distribution to Canadian retail investors.")
    print("Conclusion: This distribution is NOT compliant.")

    print("\n--- Analysis of Option C ---")
    print("Issuer: Private Issuer")
    print(f"Investor Details: Salary = ${investor_c_salary}, Net Assets = ${investor_c_net_assets}.")
    print("Exemption Analysis:")
    print("1. Private Issuer Exemption: Does not apply as the investor has 'no connection' to the company.")
    print("2. Accredited Investor Exemption: Checking financial status...")
    print(f"   - Income Test: ${investor_c_salary} is less than the ${ai_income_threshold} threshold.")
    print(f"   - Net Financial Asset Test: The problem does not provide net financial assets, but total net assets of ${investor_c_net_assets} are less than the ${ai_net_financial_assets_threshold} threshold for net *financial* assets.")
    print(f"The investor is not an accredited investor. Result: {is_accredited}.")
    print("Conclusion: This distribution is NOT compliant.")

    print("\n--- Analysis of Option D ---")
    print("Issuer: Fairstone Bank of Canada")
    print("Status: In January 2020, this entity was Fairstone Financial. It did not become a Schedule I bank until 2021.")
    print("Exemption: As it was not a bank at the time, it could not rely on the 'Financial Institution Exemption'.")
    print("Conclusion: This distribution is NOT compliant.")

    print("\n--- Analysis of Option E ---")
    print("Issuer: Caisse populaire acadienne ltée")
    print("Status: A 'caisse populaire' is a type of credit union.")
    print("Exemption: Under NI 45-106, s. 2.35, credit unions and caisses populaires can issue their own securities (like shares) without a prospectus.")
    print("The financial status of the investor (e.g., one having zero dollars of net assets) is irrelevant because this is an issuer-based exemption, not a purchaser-based one.")
    print("Conclusion: This distribution is compliant.")
    
    print("="*80)
    print("\nFinal Determination:")
    print("Both A and E appear compliant based on the Financial Institution Exemption.")
    print("However, Option E represents a textbook case for the exemption's purpose, allowing community-based cooperative financial institutions like credit unions to raise capital from their members, who are the retail public.")
    print("This scenario, including the detail about an investor with no assets, directly tests the understanding that this exemption is based on the issuer's regulated status, not the purchaser's wealth. It is the most unequivocally correct answer presented.")

if __name__ == "__main__":
    securities_regulation_analysis()