import sys
# This script is for illustrative purposes. It reflects the regulations as of January 2020.
# For any investment decision, consult with a qualified professional.

def check_accredited_investor_status():
    """
    Evaluates the investor's financial details from Option C against the
    'Accredited Investor' thresholds in Ontario's securities regulations
    (National Instrument 45-106).
    """

    # --- Financial Details from Option C ---
    investor_salary = 35000
    investor_net_assets = 10000
    
    # Financial assets are a component of net assets. In the best-case scenario for the investor,
    # we can assume all their net assets are also financial assets.
    investor_financial_assets = 10000

    # --- Accredited Investor Thresholds (for an individual) ---
    # 1. Income Test: Net income before taxes exceeds this in the 2 most recent years.
    income_test_threshold = 200000

    # 2. Financial Asset Test: Net financial assets (alone or with a spouse) exceed this.
    financial_asset_test_threshold = 1000000

    # 3. Net Asset Test: Net assets (alone or with a spouse) exceed this.
    net_asset_test_threshold = 5000000
    
    # --- Evaluation ---
    is_accredited_by_income = investor_salary >= income_test_threshold
    is_accredited_by_financial_assets = investor_financial_assets >= financial_asset_test_threshold
    is_accredited_by_net_assets = investor_net_assets >= net_asset_test_threshold

    print("Analysis of Investor in Option C against Accredited Investor Rules:")
    print("="*65)
    
    # Print Income Test Evaluation
    print("1. Income Test Evaluation:")
    print(f"   Investor's Income:      {investor_salary:10,d}")
    print(f"   Required Threshold:       {income_test_threshold:10,d}")
    print(f"   --> Meets Income Test:    {is_accredited_by_income}")
    print("-" * 65)

    # Print Financial Asset Test Evaluation
    print("2. Financial Asset Test Evaluation:")
    print(f"   Investor's Fin. Assets: {investor_financial_assets:10,d} (best case assumption)")
    print(f"   Required Threshold:       {financial_asset_test_threshold:10,d}")
    print(f"   --> Meets Fin. Asset Test: {is_accredited_by_financial_assets}")
    print("-" * 65)
    
    # Print Net Asset Test Evaluation
    print("3. Net Asset Test Evaluation:")
    print(f"   Investor's Net Assets:    {investor_net_assets:10,d}")
    print(f"   Required Threshold:       {net_asset_test_threshold:10,d}")
    print(f"   --> Meets Net Asset Test:   {is_accredited_by_net_assets}")
    print("="*65)

    if not (is_accredited_by_income or is_accredited_by_financial_assets or is_accredited_by_net_assets):
        print("\nConclusion: The investor in Option C does not qualify as an Accredited Investor.")
        print("As the investor also has no connection to the company, this distribution is non-compliant.")
    else:
        # This branch is not expected to be hit with the given numbers
        print("\nConclusion: The investor qualifies as an Accredited Investor.")

if __name__ == '__main__':
    # Execute the analysis
    check_accredited_investor_status()
    # Provide the final reasoning based on all options
    print("\nFinal Analysis of All Options:")
    print("Options A, B, and D are incorrect as general distributions of shares or bonds to the retail public by banks typically require a prospectus.")
    print("Option C is incorrect, as demonstrated by the calculations above; the investor is not accredited.")
    print("Option E is correct because the Ontario Securities Act provides a specific prospectus exemption for securities issued by a caisse populaire to its members or prospective members, regardless of their financial status.")