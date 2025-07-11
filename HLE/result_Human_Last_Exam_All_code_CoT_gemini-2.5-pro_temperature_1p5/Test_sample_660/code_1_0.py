import collections

# Define the parameters of the model.
# As per the problem description, t, p, and F are decreasing functions.
# We define them using dictionaries for the discrete values of theta and x.
t = {0: 100, 1: 50} # Tax liability is lower for truthful firms
p = {0: 0.8, 1: 0.5} # RS audit probability is lower for refusals from strict auditors
F = {0: 1000, 1: 100} # Penalty is lower for truthful firms

def cost(theta, x, t, p, F):
    """
    Calculates the expected cost C(theta, x).
    C(theta, x) = t(theta) + x * (1 - theta) * p(x) * F(theta)
    """
    # For a lenient auditor (x=0), the probability of refusal is 0.
    if x == 0:
        return t[theta]
    # For a strict auditor (x=1)
    else:
        # For a truthful firm (theta=1), the term (1-theta) is 0, so no penalty.
        if theta == 1:
            return t[theta]
        # For a malpractice firm (theta=0), there is an expected penalty.
        else:
            return t[theta] + p[x] * F[theta]

def solve():
    """
    Finds the set of theta for which firms are 'more truthful' AND 'choose lenient auditors'.
    """
    final_set = set()
    theta_values = [0, 1]

    print("Analyzing which firm types meet the criteria...\n")

    for theta in theta_values:
        # Condition 1: The company is 'more truthful'.
        is_more_truthful = (theta == 1)

        # Condition 2: The company 'chooses more lenient auditors'.
        # This means the cost with a lenient auditor (x=0) must be less than or equal to
        # the cost with a strict auditor (x=1).
        cost_lenient = cost(theta, 0, t, p, F)
        cost_strict = cost(theta, 1, t, p, F)
        chooses_lenient = (cost_lenient <= cost_strict)

        # Print the detailed analysis for each theta type
        print(f"--- Checking for theta = {theta} ---")
        if theta == 1:
            print("This firm type is 'more truthful'.")
        else:
            print("This firm type is not 'more truthful'.")

        print("Calculating costs to see if it chooses a lenient auditor:")
        # Show equation for lenient auditor
        print(f"Cost with lenient auditor (x=0): C({theta}, 0) = t({theta}) = {t[theta]}")

        # Show equation for strict auditor
        cost_strict_penalty_term = f"{1-theta} * p(1) * F({theta})"
        cost_strict_penalty_value = (1-theta) * p[1] * F[theta]
        print(f"Cost with strict auditor (x=1): C({theta}, 1) = t({theta}) + (1-{theta}) * p(1) * F({theta})")
        print(f" = {t[theta]} + {1-theta} * {p[1]} * {F[theta]} = {t[theta] + cost_strict_penalty_value}")

        if chooses_lenient:
            print(f"Result: Since {cost_lenient} <= {cost_strict}, the firm will choose the lenient auditor.")
        else:
             print(f"Result: Since {cost_lenient} > {cost_strict}, the firm will not choose the lenient auditor.")

        # Check if both conditions are met
        if is_more_truthful and chooses_lenient:
            final_set.add(theta)
            print("Conclusion: Theta={theta} satisfies BOTH conditions.\n".format(theta=theta))
        else:
            print("Conclusion: Theta={theta} does not satisfy both conditions.\n".format(theta=theta))

    print("-----------------------------------------")
    print(f"The final set of values for theta is: {final_set}")

solve()

<<<list(set([1]))[0]>>>