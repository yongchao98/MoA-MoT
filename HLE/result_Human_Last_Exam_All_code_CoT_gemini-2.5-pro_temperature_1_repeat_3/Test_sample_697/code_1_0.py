def check_accredited_investor_status(salary, net_assets):
    """
    Checks if an individual investor meets the key financial tests
    to be considered an 'Accredited Investor' in Ontario.
    This is a simplified check of the most common criteria.
    """
    
    # Define the financial thresholds for an individual Accredited Investor
    income_threshold = 200000
    net_asset_threshold = 5000000
    financial_asset_threshold = 1000000

    print("Analyzing investor status based on provided financials from Option C...")
    print(f"Investor's Salary: ${salary:,.2f}")
    print(f"Investor's Net Assets: ${net_assets:,.2f}")
    print("-" * 30)
    print("Checking against Accredited Investor thresholds:")
    
    # Check 1: Income Test
    is_income_qualified = salary > income_threshold
    print(f"1. Income Test: Is salary of ${salary:,.2f} > ${income_threshold:,.2f}? -> {is_income_qualified}")
    
    # Check 2: Net Asset Test
    is_net_asset_qualified = net_assets > net_asset_threshold
    print(f"2. Net Asset Test: Is net assets of ${net_assets:,.2f} > ${net_asset_threshold:,.2f}? -> {is_net_asset_qualified}")
    
    # Check 3: Financial Asset Test (Note: financial assets cannot exceed net assets)
    # We assume financial assets are less than or equal to net assets.
    is_financial_asset_qualified = net_assets > financial_asset_threshold
    print(f"3. Financial Asset Test: Are financial assets (assumed <= ${net_assets:,.2f}) > ${financial_asset_threshold:,.2f}? -> {is_financial_asset_qualified}")
    print("-" * 30)

    if not is_income_qualified and not is_net_asset_qualified and not is_financial_asset_qualified:
        print("Conclusion: The investor does not meet the financial requirements to be an Accredited Investor.")
        print("Therefore, the distribution in Option C without a prospectus is not compliant.")
    else:
        print("Conclusion: The investor meets at least one financial requirement to be an Accredited Investor.")

# Financial details from Option C
investor_salary = 35000
investor_net_assets = 10000

# Run the check
check_accredited_investor_status(investor_salary, investor_net_assets)
