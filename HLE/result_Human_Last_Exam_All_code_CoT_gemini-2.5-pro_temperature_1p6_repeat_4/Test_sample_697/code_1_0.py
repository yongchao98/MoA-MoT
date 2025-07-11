# Investor's financial details from Option C
investor_salary = 35000
investor_net_assets = 10000
# Financial assets are not given, but they cannot exceed net assets.
# We can assume investor_financial_assets <= 10000 for this check.
investor_financial_assets = 10000 

# Accredited Investor thresholds for an individual in Ontario (as of Jan 2020)
income_threshold = 200000
financial_assets_threshold = 1000000
net_assets_threshold = 5000000

# --- Check Compliance with Accredited Investor Exemption ---

# Check 1: Income Test
is_income_met = investor_salary > income_threshold
# Check 2: Financial Asset Test
is_financial_assets_met = investor_financial_assets > financial_assets_threshold
# Check 3: Net Asset Test
is_net_assets_met = investor_net_assets >= net_assets_threshold

is_accredited = is_income_met or is_financial_assets_met or is_net_assets_met

print("--- Analysis of Investor in Option C ---")
print(f"To qualify as an accredited investor, one of the following must be true:")
print(f"1. Annual Income > ${income_threshold:,}")
print(f"2. Net Financial Assets > ${financial_assets_threshold:,}")
print(f"3. Net Assets >= ${net_assets_threshold:,}")
print("\n--- Investor's Financials vs. Thresholds ---")
print(f"Income Check: Investor's salary of ${investor_salary:,} is NOT greater than the threshold of ${income_threshold:,}. Met: {is_income_met}")
print(f"Financial Asset Check: Investor's financial assets of <=${investor_financial_assets:,} is NOT greater than the threshold of ${financial_assets_threshold:,}. Met: {is_financial_assets_met}")
print(f"Net Asset Check: Investor's net assets of ${investor_net_assets:,} is NOT greater than or equal to the threshold of ${net_assets_threshold:,}. Met: {is_net_assets_met}")

print("\n--- Conclusion for Option C ---")
if not is_accredited:
    print("The investor does not qualify as an accredited investor. As no other exemption applies, this distribution is NOT compliant.")
else:
    print("The investor qualifies as an accredited investor.")
