import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_securities_distributions():
    """
    Analyzes several securities distribution scenarios based on Ontario
    securities regulations as of January 2020.
    """

    print("Analyzing compliance of securities distributions in Ontario...\n")

    # --- Analysis for Option C: Accredited Investor Test ---
    print("--- Analysis of Option C ---")
    investor_c_salary = 35000
    investor_c_net_assets = 10000

    # Financial thresholds for an individual to be an "Accredited Investor"
    income_threshold = 200000
    net_financial_assets_threshold = 1000000
    net_assets_threshold = 5000000

    print(f"Investor C has a salary of ${investor_c_salary} and net assets of ${investor_c_net_assets}.")
    
    is_accredited = False
    if investor_c_salary > income_threshold:
        is_accredited = True
    # Net financial assets are not given, but cannot exceed net assets.
    if investor_c_net_assets > net_assets_threshold:
        is_accredited = True

    print(f"Is the investor an 'Accredited Investor'? {is_accredited}.")
    print("Conclusion for C: The investor does not meet the financial thresholds (e.g., income > $200k or net assets > $5M). Since there is no other connection to the company, the distribution is NON-COMPLIANT.\n")

    # --- Explanatory Analysis for Other Options ---
    print("--- Analysis of Other Options ---")

    print("Option A (Bank of China Bonds): COMPLIANT. A prospectus exemption (NI 45-106, s. 2.38(a)) exists for the distribution of bonds by a regulated Canadian bank (Schedule I, II, or III).\n")
    
    print("Option B (JPMorgan Chase Bonds): NON-COMPLIANT. JPMorgan Chase is not a bank under the Canadian Bank Act, so the specific bank exemption does not apply. A distribution to a large number of Canadian investors requires a prospectus.\n")
    
    print("Option D (Fairstone Bank Shares): NON-COMPLIANT. The prospectus exemption for Canadian banks (NI 45-106, s. 2.38(a)) applies to debt and deposit instruments (like bonds), but NOT to shares.\n")
    
    print("Option E (Caisse Populaire Shares): COMPLIANT. A prospectus exemption (NI 45-106, s. 2.38(b)) exists for distributions of securities by a caisse populaire to its members. The compliance is based on the investor's membership status, not their financial wealth. The fact that a member has zero assets does not invalidate the use of this exemption.\n")

    print("--- Final Conclusion ---")
    print("Based on the analysis, Option E represents a clear and valid use of a prospectus exemption under Ontario securities law.")


# Run the analysis
analyze_securities_distributions()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())