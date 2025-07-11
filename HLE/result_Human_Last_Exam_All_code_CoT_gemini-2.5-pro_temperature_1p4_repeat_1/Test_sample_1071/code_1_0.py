def solve_bounds():
    """
    This function calculates the lower and upper bounds of t based on the
    mathematical derivation explained above.
    """

    # The problem reduces to finding t such that for R = [-2, 2t], f(R) is a subset of R, where f(x) = 1/x.
    # This leads to two inequalities that determine the bounds of t.

    # 1. Inequality from comparing the upper bounds: -1/2 <= 2*t
    #    This sets the lower bound for t.
    n1 = -1.0
    d1 = 2.0
    c1 = 2.0
    lower_bound = (n1 / d1) / c1
    print("Derivation of the lower bound for t:")
    print(f"The inequality is {n1}/{d1} <= {c1}*t.")
    print(f"Solving for t gives t >= {lower_bound}.")
    print("-" * 20)

    # 2. Inequality from comparing the lower bounds: 1/(2t) >= -2
    #    This can be rewritten as 1 <= -4*t (for t < 0).
    #    This sets the upper bound for t.
    n2 = 1.0
    c2 = -4.0
    upper_bound = n2 / c2
    print("Derivation of the upper bound for t:")
    print("The inequality is 1/(2t) >= -2, which simplifies to 1 <= -4t.")
    print(f"Solving for t gives t <= {n2}/{c2} = {upper_bound}.")
    print("-" * 20)
    
    print(f"Since t must be both >= {lower_bound} and <= {upper_bound}, the only possible value is t = -0.25.")
    print(f"Final Answer: The lower bound for t is {lower_bound} and the upper bound is {upper_bound}.")

solve_bounds()