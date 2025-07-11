import sys
import io

# Redirect stdout to capture print statements
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def analyze_securities_distributions():
    """
    Analyzes five scenarios against Ontario securities regulations (as of Jan 2020)
    to determine which one is compliant.
    """

    # Scenario C financial details
    c_salary = 35000
    c_net_assets = 10000

    # Scenario D/E investor details
    d_e_investor_net_assets = 0
    d_e_investor_is_unemployed = True
    
    # Accredited Investor Thresholds (simplified)
    AI_INCOME_THRESHOLD = 200000
    AI_NET_ASSETS_THRESHOLD = 5000000
    AI_NET_FINANCIAL_ASSETS_THRESHOLD = 1000000

    print("Analyzing compliance of securities distributions in Ontario (as of Jan 2020)...")
    print("-" * 70)

    # --- Scenario A ---
    is_a_compliant = False
    reason_a = ("A: Bank of China (Canada) (Schedule III Bank) distributing bonds (debt) without a prospectus.\n"
                "   - Rule Check: The exemption for banks (NI 45-106 s. 2.35) applies to debt securities from Schedule I, II, or III banks.\n"
                "   - Analysis: This fits the rule. However, retail fundraising is unusual for a Schedule III bank and might conflict with banking regulations, making its overall compliance status less certain than other options.")
    print(reason_a)
    print("-" * 70)

    # --- Scenario B ---
    is_b_compliant = False
    reason_b = ("B: JPMorgan Chase (not a Canadian Schedule I/II/III bank) distributing bonds without a prospectus.\n"
                "   - Rule Check: The bank exemption only applies to banks chartered under the Bank Act (Canada).\n"
                "   - Analysis: Not compliant, as JPMorgan Chase is not a Schedule I, II, or III bank.")
    print(reason_b)
    print("-" * 70)

    # --- Scenario C ---
    is_c_compliant = False
    is_accredited_investor = (c_salary > AI_INCOME_THRESHOLD or c_net_assets > AI_NET_ASSETS_THRESHOLD)
    reason_c = (f"C: Private issuer distributing shares to an unconnected investor with a salary of ${c_salary:,} and net assets of ${c_net_assets:,}.\n"
                "   - Rule Check: Private issuer exemption requires a close connection. The accredited investor exemption requires meeting financial thresholds.\n"
                "   - Analysis: The investor has no connection and does not meet the accredited investor thresholds (e.g., income must be > ${AI_INCOME_THRESHOLD:,}). Not compliant.")
    print(reason_c)
    print("-" * 70)

    # --- Scenario D ---
    is_d_compliant = False
    reason_d = ("D: Fairstone Bank of Canada (Schedule I Bank) distributing shares (equity) without a prospectus.\n"
                "   - Rule Check: The exemption for banks (NI 45-106 s. 2.35) applies to DEBT, not EQUITY.\n"
                "   - Analysis: Not compliant, as the wrong type of security is being issued under this exemption.")
    print(reason_d)
    print("-" * 70)

    # --- Scenario E ---
    is_e_compliant = True
    reason_e = (f"E: Caisse populaire acadienne ltée distributing shares without a prospectus to retail investors, one with ${d_e_investor_net_assets} of net assets.\n"
                "   - Rule Check: A specific exemption (NI 45-106 s. 2.36) exists for credit unions and caisses populaires to issue their own securities, including shares, to their community.\n"
                "   - Analysis: This is compliant. The exemption is designed for this exact purpose, and the financial status of individual retail investors is irrelevant to its application.")
    print(reason_e)
    print("-" * 70)

    final_answer = 'E'
    print(f"Final Conclusion: Scenario {final_answer} describes the most clearly and robustly compliant distribution.")


# Run the analysis
analyze_securities_distributions()

# Restore stdout and print the captured output to the actual console
sys.stdout = old_stdout
print(captured_output.getvalue())

# Final Answer Marker
print("<<<E>>>")