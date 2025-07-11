def solve_puzzle():
    """
    Calculates the maximal probability p for Alice's success.
    """
    # Total number of boxes
    total_boxes = 20

    # Alice can open k boxes, where 1 <= k < total_boxes.
    # To maximize her win probability, Alice should open the maximum number of boxes possible.
    k = total_boxes - 1

    # The winning probability is given by the formula p(k) = 1 - 2 / (k + 1).
    # We calculate this for k = 19.
    
    numerator = 2
    denominator = k + 1
    
    # We want to print the equation step-by-step
    print(f"Alice's optimal strategy is to open the maximum number of boxes allowed, k = {k}.")
    print("The probability of success p is calculated with the formula: p = 1 - 2 / (k + 1)")
    print("Substituting k with its value:")
    print(f"p = 1 - {numerator} / ({k} + 1)")
    print(f"p = 1 - {numerator} / {denominator}")
    print(f"p = 1 - 1 / {denominator // numerator}")
    print(f"p = {denominator // numerator - 1} / {denominator // numerator}")

solve_puzzle()