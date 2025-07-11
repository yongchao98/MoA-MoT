def find_root_cause_for_rmbs_value():
    """
    This script analyzes the factors determining non-agency RMBS value
    during the 2004-2008 US housing crisis to identify the root cause.
    """

    print("Analyzing the root cause factor for non-agency RMBS value (2004-2008)...")
    print("-" * 70)
    
    # Step 1: Identify factors that are outcomes or symptoms, not root causes.
    print("Evaluating Factor E (Default rates) and G (Recovery rates):")
    print("These are direct inputs into a security's value, but they are the *result* of risk, not the origin of the risk. They are symptoms of the problem.")
    print("\nEvaluating Factor H (The rating of the credit agency):")
    print("The ratings from agencies proved to be inaccurate and unreliable. The market's realization of this fact is what caused the crash. Therefore, the rating was a flawed signal, not a determinant of true value.")
    
    print("-" * 70)
    
    # Step 2: Evaluate underlying risk characteristics.
    print("Evaluating Factor C (Average FICO scores):")
    print("Low FICO scores were a clear indicator of high credit risk in the mortgage pools. This is a very strong factor.")
    
    print("-" * 70)

    # Step 3: Identify the most fundamental cause that explains the other factors.
    print("Evaluating Factor F (The quality of the loan issuer and RMBS originator):")
    print("This is the most fundamental 'root cause'. The originator's underwriting standards (or lack thereof) are the reason the mortgage pools contained high-risk elements in the first place.")
    print("The decision to approve loans for low-FICO borrowers (Factor C) and to create other risky loan products was made by the originator.")
    print("Therefore, the originator's quality is the source of the risk that ultimately determined the security's value.")

    print("-" * 70)
    print("Conclusion: The quality and practices of the loan issuer and originator are the most upstream factor and the primary root cause.")
    
    # Final Answer Output
    print("\nFinal Answer Selection:")
    print("F. The quality of the loan issuer and RMBS originator")
    print("<<<F>>>")

find_root_cause_for_rmbs_value()