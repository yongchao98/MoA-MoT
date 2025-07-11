# This is a conceptual model to illustrate the causal relationships
# determining non-agency RMBS value during the 2004-2008 crisis.
# It represents the logic rather than performing a real financial calculation.

def analyze_rmbs_root_cause():
    """
    Analyzes the factors contributing to the RMBS crisis and identifies the root cause.
    """
    # Define the factors from the answer choices
    factors = {
        "A": "Stock market level of the S&P500 (Consequence/Correlated)",
        "B": "The percent of floating rate debt in the pool (Loan Characteristic)",
        "C": "Average FICO scores on the loans (Loan Characteristic)",
        "D": "The 10 year US Treasury rate (Market Benchmark)",
        "E": "Default rates (Direct Cause/Outcome)",
        "F": "The quality of the loan issuer and RMBS originator (Root Cause)",
        "G": "Recovery rates (Direct Cause/Outcome)",
        "H": "The rating of the credit agency on issuance (Enabling Factor)"
    }

    # The causal chain starts with the originator.
    # In the 2004-2008 period, originator quality was the key issue.
    root_cause_factor_quality = "Low"
    root_cause_factor_key = "F"

    print("Step 1: Identify the Root Cause")
    print(f"The analysis begins with the most fundamental factor: {factors[root_cause_factor_key]}.")
    print(f"During the 2004-2008 period, the quality of many originators was '{root_cause_factor_quality}' due to perverse incentives.\n")


    # Low-quality origination leads to risky loan characteristics.
    print("Step 2: Trace the Impact of the Root Cause")
    if root_cause_factor_quality == "Low":
        loan_characteristic_fico = "Low"  # Factor C
        loan_characteristic_floating_rate = "High" # Factor B
        print(f"Low-quality origination practices led to creating pools of loans with risky characteristics:")
        print(f" - Low Average FICO scores (Factor C)")
        print(f" - High percentage of floating-rate debt (Factor B)\n")

    # Risky loan characteristics lead to poor outcomes.
    print("Step 3: Observe the Final Outcomes")
    if loan_characteristic_fico == "Low" and loan_characteristic_floating_rate == "High":
        outcome_defaults = "High" # Factor E
        outcome_recoveries = "Low"  # Factor G
        print(f"These risky loan pools inevitably resulted in poor performance:")
        print(f" - High Default Rates (Factor E)")
        print(f" - Low Recovery Rates (Factor G) as mass foreclosures crashed the housing market.\n")


    # Final conclusion based on the causal chain.
    print("Final Conclusion Equation:")
    final_equation = f"({factors[root_cause_factor_key]}) -> (Low {factors['C']} & High {factors['B']}) -> (High {factors['E']} & Low {factors['G']}) -> RMBS Value Collapse"
    print(f"The logical flow shows the chain of causation: {final_equation}")
    print("\nTherefore, the quality of the loan issuer and RMBS originator was the root cause factor most likely to determine the value of non-agency RMBS.")

# Execute the analysis
analyze_rmbs_root_cause()