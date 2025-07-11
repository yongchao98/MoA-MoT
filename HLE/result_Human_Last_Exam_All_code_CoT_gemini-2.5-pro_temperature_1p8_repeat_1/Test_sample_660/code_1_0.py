def solve():
    """
    This function derives and prints the condition under which more truthful
    companies might choose more lenient auditors.

    Let C(theta, x) be the cost for a firm of type theta hiring an auditor of type x.
    theta=1 is a truthful firm, theta=0 is a firm with malpractice.
    x=1 is a strict auditor, x=0 is a lenient auditor.

    The cost is given by C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta).

    For a lenient auditor (x=0): C(theta, 0) = t(theta).
    For a strict auditor (x=1): C(theta, 1) = t(theta) + (1 - theta) * p(1) * F(theta).

    A firm's preference for a lenient auditor is measured by the cost difference
    L(theta) = C(theta, 1) - C(theta, 0).
    L(theta) = (1 - theta) * p(1) * F(theta).

    For a truthful firm (theta=1): L(1) = (1 - 1) * p(1) * F(1) = 0.
    For a malpractice firm (theta=0): L(0) = (1 - 0) * p(1) * F(0) = p(1) * F(0).

    The condition "more truthful firms choose more lenient auditors" implies
    that the preference for leniency is higher for theta=1 than for theta=0.
    Mathematically, this means L(1) >= L(0).

    Substituting the values, we get:
    0 >= p(1) * F(0)

    Since probability p(1) and penalty F(0) must be non-negative, the only way
    to satisfy this condition is for the product to be zero.
    """
    
    # Define the variables from the problem for the equation
    p_1 = "p(1)" # Probability of RS audit if strict auditor (x=1) refuses
    F_0 = "F(0)" # Penalty for a malpractice firm (theta=0)
    
    # The final equation representing the condition
    print("The condition for which companies keeping more truthful accounts may choose more lenient auditors is:")
    # Using 'times' to represent multiplication for clarity
    final_equation = f"{p_1} * {F_0} = 0"
    
    # We now print each number and symbol in the final equation.
    print(f"p(1) * F(0) = 0")


solve()
