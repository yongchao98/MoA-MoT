def analyze_rmbs_factors():
    """
    Analyzes the factors affecting non-agency RMBS value from 2004-2008
    to determine the root cause.
    """

    print("Step 1: Understanding the difference between a proximate cause and a root cause.")
    print("A proximate cause is the most immediate reason for an event (e.g., high defaults).")
    print("A root cause is the fundamental reason that created the conditions for the proximate cause to occur.")
    print("-" * 70)

    print("Step 2: Evaluating the answer choices.")
    print(" - E. Default rates and G. Recovery rates are PROXIMATE causes. They directly determine the cash flow and losses from the mortgage pool, but they don't explain *why* defaults and losses were so high.")
    print(" - H. The credit agency ratings were a major failure of oversight that enabled the crisis by misrepresenting risk, but they were not the source of the risk itself.")
    print(" - B. Floating rate debt and C. Average FICO scores are characteristics of the underlying loans. They are symptoms of poor loan quality, not the overarching cause.")
    print(" - F. The quality of the loan issuer and RMBS originator gets to the core of the problem. The 'originate-to-distribute' model incentivized originators to create and sell as many loans as possible with no regard for quality (e.g., 'NINJA' loans).")
    print("-" * 70)

    print("Step 3: Identifying the root cause.")
    print("The systemic breakdown in lending standards by loan originators is the fundamental or 'root' cause. This poor origination quality is the reason that default rates skyrocketed once housing prices stalled.")
    
    # In this qualitative problem, there is no numerical equation.
    # We are identifying the letter corresponding to the root cause factor.
    final_answer = "F"
    
    # The final output presents the conclusion of the analysis.
    print(f"\nConclusion: The root cause factor is '{final_answer}'. This corresponds to 'The quality of the loan issuer and RMBS originator'.")

    # The final answer is printed in the required format.
    # The instruction to "output each number in the final equation" is not applicable
    # to this qualitative problem, so we present the concluding letter.
    print(f"<<<{final_answer}>>>")

# Execute the analysis
analyze_rmbs_factors()