def solve_auditor_choice_problem():
    """
    This function analyzes the auditor choice problem and prints the derivation
    for the condition under which more truthful firms choose more lenient auditors.
    """
    # Represent the parameters symbolically using strings
    t_0, t_1 = "t(0)", "t(1)"
    p_1 = "p(1)"
    F_0 = "F(0)"

    # Define the cost for the malpractice firm (theta=0) for each auditor choice
    cost_malpractice_lenient = t_0
    cost_malpractice_strict = f"{t_0} + {p_1} * {F_0}"

    # Define the cost for the truthful firm (theta=1) for each auditor choice
    cost_truthful_lenient = t_1
    cost_truthful_strict = t_1 # The penalty term becomes 1*(1-1)*p(1)*F(1) = 0

    # Print the step-by-step reasoning
    print("Step 1: Define the cost functions for each firm type and auditor choice.")
    print("The cost C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta).")
    print("-" * 50)
    print("For a malpractice firm (theta=0):")
    print(f"  Cost with lenient auditor (x=0): C(0,0) = {cost_malpractice_lenient}")
    print(f"  Cost with strict auditor (x=1):   C(0,1) = {cost_malpractice_strict}")
    print("\nFor a truthful firm (theta=1):")
    print(f"  Cost with lenient auditor (x=0): C(1,0) = {cost_truthful_lenient}")
    print(f"  Cost with strict auditor (x=1):   C(1,1) = {cost_truthful_strict}")
    print("-" * 50)
    
    print("Step 2: Determine the conditions for the required choices.")
    print("We need the truthful firm (theta=1) to choose the lenient auditor (x=0) and")
    print("the malpractice firm (theta=0) to choose the strict auditor (x=1).")
    
    print("\nCondition for the truthful firm to choose lenient (C(1,0) <= C(1,1)):")
    print(f"  {cost_truthful_lenient} <= {cost_truthful_strict}, which simplifies to {t_1} <= {t_1}.")
    print("  This is always true. The truthful firm is always indifferent.")
    
    print("\nCondition for the malpractice firm to choose strict (C(0,1) <= C(0,0)):")
    print(f"  {cost_malpractice_strict} <= {cost_malpractice_lenient}")
    print(f"  This simplifies to {p_1} * {F_0} <= 0.")
    print("-" * 50)
    
    print("Step 3: Conclude the final condition.")
    print("Since p(1) (a probability) and F(0) (a penalty) must be non-negative,")
    print(f"the inequality '{p_1} * {F_0} <= 0' can only hold if they are equal to zero.")
    
    print("\nTherefore, the condition under which it is possible for more truthful firms")
    print("to choose more lenient auditors than less truthful firms is:")
    
    # Print each symbol in the final equation as requested
    print(p_1, "*", F_0, "=", "0")

# Run the analysis
solve_auditor_choice_problem()