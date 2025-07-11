def solve_marble_problem():
    """
    Calculates the probability that the marble escapes.

    The problem can be solved by considering the ratio of escape vs. melt probabilities.
    For a symmetric random walk starting at n with absorbing states at a and b,
    the ratio of probabilities of escaping (hitting b first) to melting (hitting a first)
    is given by:
    P(escape) / P(melt) = 3^(|a-n| - |b-n|)

    Let p be the probability of escape. Then 1-p is the probability of melting.
    p / (1-p) = 3^(|a-n| - |b-n|)
    Solving for p gives:
    p = 1 / (1 + 3^(|b-n| - |a-n|))
    """
    n = 0      # Starting bin
    a = 2024   # Torch bin (melt)
    b = 2025   # Portal bin (escape)

    # Calculate the exponent
    exp_val = abs(b - n) - abs(a - n)

    # Calculate the final probability
    probability = 1 / (1 + 3**exp_val)

    print("Let p be the probability of the marble escaping.")
    print("The ratio of escape to melt probability is determined by the formula:")
    print(f"p / (1 - p) = 3^(|{a} - {n}| - |{b} - {n}|)")
    print("\nSolving for p, we get:")
    print(f"p = 1 / (1 + 3^(|{b} - {n}| - |{a} - {n}|))")
    
    print("\nPlugging in the values:")
    print(f"The exponent is |{b} - {n}| - |{a} - {n}| = {abs(b-n)} - {abs(a-n)} = {exp_val}")
    print(f"So, p = 1 / (1 + 3^{exp_val})")
    print(f"p = 1 / (1 + {3**exp_val})")
    print(f"p = 1 / {1 + 3**exp_val}")
    print(f"p = {probability}")

solve_marble_problem()