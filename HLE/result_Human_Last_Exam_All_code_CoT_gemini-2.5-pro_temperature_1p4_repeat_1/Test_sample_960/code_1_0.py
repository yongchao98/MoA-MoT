def solve():
    """
    Calculates the probability that the number of good and defective products will ever be equal.
    """
    # Initial number of good (white) and defective (black) products
    w = 2
    b = 1

    print(f"Initial state: W_0 = {w}, B_0 = {b}")

    # The problem asks for the probability of reaching a state where the number of good and defective products are equal.
    # This is a classic problem in the theory of Polya's Urns.
    # Let P(w, b) be the probability of reaching equality, starting with w good and b defective items.
    # We assume w >= b without loss of generality.
    if w < b:
        w, b = b, w  # Ensure w is the larger count

    print(f"We are looking for the probability of reaching equality, starting from w={w} and b={b}.")
    
    # The formula for the probability of reaching equality for a Polya's Urn starting with w > b is:
    # P = 1 if w <= b + 1
    # P = b / (w - 1) if w > b + 1

    print("Checking the condition: w <= b + 1")
    print(f"{w} <= {b} + 1")
    print(f"{w} <= {b + 1}")
    
    if w <= b + 1:
        probability = 1.0
        print("The condition is true.")
        print("Therefore, the probability is 1.")
    else:
        probability = b / (w - 1)
        print("The condition is false.")
        print(f"The probability is calculated as b / (w - 1) = {b} / ({w} - 1) = {probability}")

    # The question asks for the upper bound. Since we found the exact probability,
    # this value is the least upper bound.
    print(f"\nThe exact probability of reaching a 50% good product ratio is {probability}.")
    print(f"Thus, the least upper bound for this probability is {probability}.")
    final_equation = f"{w} <= {b} + 1 -> True, Probability = 1.0"
    print(f"\nFinal calculation summary: {final_equation}")


solve()