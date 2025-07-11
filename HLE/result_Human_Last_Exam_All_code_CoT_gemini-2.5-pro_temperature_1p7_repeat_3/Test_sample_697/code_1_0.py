import sys

# Define constants for Accredited Investor tests (as of Jan 2020)
AI_NET_INCOME_TEST = 200000
AI_NET_ASSET_TEST = 5000000
AI_NET_FINANCIAL_ASSET_TEST = 1000000

# Data for the investor in Option C
investor_c_salary = 35000
investor_c_net_assets = 10000

def is_accredited_investor(salary, net_assets):
    """
    Checks if an investor meets the simplified financial tests to be an 'accredited investor'.
    """
    if salary > AI_NET_INCOME_TEST:
        return True
    if net_assets >= AI_NET_ASSET_TEST:
        return True
    # Assuming net_assets could be considered net_financial_assets for this test
    if net_assets > AI_NET_FINANCIAL_ASSET_TEST:
        return True
    return False

def analyze_scenarios():
    """
    Analyzes each scenario based on Ontario securities regulations.
    """
    print("Analyzing Compliance with Ontario Securities Regulations (as of Jan 2020)")
    print("="*70)

    # --- Option A ---
    print("Analysis for A: Bonds from Bank of China (Canada) [Schedule III Bank]")
    print("Rule: NI 45-106 s. 2.35 exempts securities from a bank.")
    print("Potential Issue: The exemption does not apply to subordinated bonds unless they have a specific high credit rating. The status of the bonds is not specified.")
    print("Compliance: Uncertain. If the bonds are subordinated and unrated, it would NOT be compliant.")
    print("-" * 70)

    # --- Option B ---
    print("Analysis for B: Bonds from JPMorgan Chase [US Parent Co.]")
    print("Rule: Issuer is not a bank under the Canadian Bank Act.")
    print("Potential Issue: The bank exemption does not apply. Without another applicable exemption (like 'accredited investor'), a prospectus is required for distribution to the public.")
    print("Compliance: NOT Compliant.")
    print("-" * 70)

    # --- Option C ---
    print(f"Analysis for C: Shares from Private Issuer to investor with:")
    print(f" * Salary: {investor_c_salary}")
    print(f" * Net Assets: {investor_c_net_assets}")
    investor_is_accredited = is_accredited_investor(investor_c_salary, investor_c_net_assets)
    print("Rule: Distribution from a private issuer without a prospectus is generally restricted to connected persons or accredited investors.")
    print(f"Is investor accredited based on financial details? {investor_is_accredited}.")
    print("Potential Issue: The investor has no connection and does not meet the income or asset tests to be an accredited investor.")
    print("Compliance: NOT Compliant.")
    print("-" * 70)

    # --- Option D ---
    print("Analysis for D: Shares from Fairstone Bank of Canada [Schedule I Bank]")
    print("Rule: NI 45-106 s. 2.35 provides a clear prospectus exemption for securities issued by a Schedule I bank. This applies to shares.")
    print("Note: The individual investor's financial status (e.g., zero net assets) is a matter for a dealer's 'suitability' obligation, but it does not invalidate the prospectus exemption for the distribution itself.")
    print("Compliance: Compliant.")
    print("-" * 70)
    
    # --- Option E ---
    print("Analysis for E: Shares from Caisse populaire acadienne ltée [Credit Union]")
    print("Rule: NI 45-106 s. 2.36 provides a prospectus exemption for securities issued by a credit union.")
    print("Note: Like option D, investor financial status is a suitability matter. The rules are sound, making this a strong candidate for compliance.")
    print("Compliance: Compliant.")
    print("=" * 70)
    
    print("\nConclusion:")
    print("Options B and C are clearly non-compliant.")
    print("Option A is questionable because of the potential for the bonds to be subordinated, which would void the exemption.")
    print("Both D and E describe distributions that are compliant under NI 45-106. However, the regulatory framework for federal banks (D) is highly uniform across Canada, making it a very clear-cut example of a compliant distribution. Therefore, it stands out as the best and most unambiguous answer.")

if __name__ == "__main__":
    analyze_scenarios()