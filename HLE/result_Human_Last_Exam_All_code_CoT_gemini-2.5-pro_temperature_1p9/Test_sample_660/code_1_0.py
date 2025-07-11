from sympy import symbols, aprint

def solve_auditor_choice():
    """
    This function derives and prints the condition on a firm's truthfulness parameter (theta)
    for choosing a lenient auditor, based on the model developed.
    """
    
    # Define symbols for the parameters in the model.
    # theta: firm's truthfulness
    # p(1): probability of government audit after refusal by a strict auditor
    # q(0): probability of government audit after being cleared by a lenient auditor
    # q(1): probability of government audit after being cleared by a strict auditor
    theta = symbols('theta')
    p_1, q_0, q_1 = symbols('p(1) q(0) q(1)')
    
    # The derivation shows that firms with theta greater than a specific threshold
    # will choose the lenient auditor. This threshold is a function of the audit probabilities.
    # The condition for choosing a lenient auditor is:
    # theta > (p(1) - q(0)) / (p(1) - q(1))
    
    # We print the final equation that defines the set of values for theta.
    # The problem implies parameter values exist such that this sorting occurs.
    # The necessary (and counter-intuitive) conditions are p(1) < q(0) < q(1).
    
    equation_text = f"theta > (p(1) - q(0)) / (p(1) - q(1))"
    
    print("The set of values of theta for which companies choose more lenient auditors is given by the inequality:")
    print(equation_text)
    
    # We can also represent this with sympy, though printing the text is clearer.
    # threshold = (p_1 - q_0) / (p_1 - q_1)
    # inequality = theta > threshold
    # print("\nSymbolic representation:")
    # aprint(inequality)


solve_auditor_choice()