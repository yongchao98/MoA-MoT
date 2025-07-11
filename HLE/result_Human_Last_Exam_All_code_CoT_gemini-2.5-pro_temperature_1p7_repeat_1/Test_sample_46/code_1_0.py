import sys

def analyze_rmbs_factors():
    """
    This function analyzes the factors determining non-agency RMBS value
    from 2004 to 2008 and prints a logical conclusion.
    """
    print("Task: Identify the root cause factor for non-agency RMBS value from 2004-2008.")
    print("-" * 70)
    print("Step 1: Understand RMBS Valuation.")
    print("The value of an RMBS is the present value of the cash flows from the underlying mortgages.")
    print("Conceptual Equation 1: RMBS Value = PV(Sum of Expected Monthly Payments from Loan Pool)")
    print("\n")

    print("Step 2: Identify what drives the 'Expected Monthly Payments'.")
    print("The payments depend on the performance of the underlying loans. The key variables are default rates and recovery rates upon default.")
    print("However, default rates and recovery rates (Choices E and G) are OUTCOMES, not root causes.")
    print("The root cause is the factor that *determines* the future rate of defaults.")
    print("\n")

    print("Step 3: Pinpoint the source of high future defaults.")
    print("High defaults are caused by poor quality loans being made in the first place.")
    print("This involves factors like low FICO scores (C) and risky loan structures (B).")
    print("\n")

    print("Step 4: Identify the root cause of poor quality loans.")
    print("The reason so many poor-quality loans were created was the breakdown of lending standards by their creators.")
    print("Therefore, the 'quality of the loan issuer and RMBS originator' (Choice F) is the true root cause.")
    print("Lax and often fraudulent origination practices (e.g., 'NINJA loans') created pools of loans that were destined to fail.")
    print("-" * 70)
    
    print("Final Causal Chain Equation:")
    print("Poor Issuer/Originator Quality (F) leads to")
    print("--> Bad Loans (Low FICO, high LTV, etc.) leads to")
    print("--> High Default Rates (E) leads to")
    print("--> Massive Loss of RMBS Value")
    
analyze_rmbs_factors()
<<<F>>>