import sys

# Disable writing __pycache__ files
sys.dont_write_bytecode = True

def analyze_securities_regulations():
    """
    Analyzes several securities distribution scenarios based on Ontario
    securities regulations as of January 2020 to find the compliant one.
    """
    print("Analyzing scenarios based on Ontario securities regulations (as of Jan 2020)...\n")

    # --- Scenario A ---
    print("--- Analyzing Scenario A ---")
    print("A distribution of bonds without a prospectus by Bank of China (Canada) aimed at raising funds from retail investors.")
    # Regulation Check: Bank of China (Canada) is a Schedule III bank (foreign bank subsidiary).
    # Prospectus exemptions for Schedule III banks are narrower than for domestic (Schedule I) banks.
    # A broad distribution of securities to retail investors generally requires a prospectus.
    is_A_compliant = False
    print("Result: NOT COMPLIANT. A prospectus is generally required for a Schedule III bank to distribute securities to retail investors.\n")

    # --- Scenario B ---
    print("--- Analyzing Scenario B ---")
    print("A distribution of bonds by JPMorgan Chase without a prospectus aimed at raising funds from a large number of investors in Canada.")
    # Regulation Check: JPMorgan Chase is a foreign bank. Unless the "large number of investors"
    # are all accredited investors (which is not stated), a public distribution requires a prospectus.
    is_B_compliant = False
    print("Result: NOT COMPLIANT. A distribution to a 'large number of investors' implies a public offering, which requires a prospectus.\n")

    # --- Scenario C ---
    print("--- Analyzing Scenario C ---")
    print("A distribution of shares by a private issuer to an individual investor who has no connection to the company and has disclosed the following financial details measured in Canadian dollars:")
    salary = 35000
    net_assets = 10000
    print(f"* A salary of {salary}")
    print(f"* Net assets of {net_assets}")
    # Regulation Check: Since the investor has no connection, we check if they are an "Accredited Investor".
    # Key Accredited Investor Thresholds (Individual):
    # 1. Net income > $200,000 annually
    # 2. Net financial assets > $1,000,000
    # 3. Net assets > $5,000,000
    income_threshold = 200000
    net_asset_threshold = 5000000
    is_accredited = (salary > income_threshold) or (net_assets > net_asset_threshold)
    is_C_compliant = is_accredited
    print(f"\nChecking Accredited Investor Status:")
    print(f"Investor's Salary: ${salary}. Is it > ${income_threshold}? {salary > income_threshold}")
    print(f"Investor's Net Assets: ${net_assets}. Is it > ${net_asset_threshold}? {net_assets > net_asset_threshold}")
    print("The investor does not meet the income or net asset tests to be an accredited investor.")
    print("Result: NOT COMPLIANT. The distribution is to a non-accredited investor with no other qualifying connection to the private issuer.\n")

    # --- Scenario D ---
    print("--- Analyzing Scenario D ---")
    print("A distribution of shares without a prospectus by Fairstone Bank of Canada aimed at raising funds from a large number of retail investors...")
    # Regulation Check: Fairstone Bank is a Schedule I bank. The prospectus exemption for banks primarily
    # applies to debt instruments (like GICs or bonds), not equity (shares) sold to the general retail public.
    is_D_compliant = False
    print("Result: NOT COMPLIANT. The prospectus exemption for banks typically does not cover distributions of equity (shares) to the retail public.\n")

    # --- Scenario E ---
    print("--- Analyzing Scenario E ---")
    print("A distribution of shares without a prospectus by Caisse populaire acadienne ltée aimed at raising funds from a large number of retail investors...")
    # Regulation Check: National Instrument 45-106 provides a specific and broad prospectus exemption for securities
    # issued by a credit union or caisse populaire. This allows them to raise capital from their members and the public.
    # The financial status of the individual investors is not a condition of this particular exemption.
    is_E_compliant = True
    print("Result: COMPLIANT. There is a specific prospectus exemption for securities issued by a caisse populaire (credit union).\n")

    # --- Final Conclusion ---
    if is_E_compliant:
        final_answer = 'E'
    else:
        # Fallback in case logic is flawed, though it shouldn't be.
        final_answer = "Error in analysis"

    print("-------------------------")
    print(f"The compliant scenario is E.")
    print(f"<<<{final_answer}>>>")


if __name__ == "__main__":
    analyze_securities_regulations()