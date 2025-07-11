def solve_snail_problem():
    """
    Calculates the maximal distance the snail could have traveled.

    The problem can be solved by considering an optimal arrangement of observers
    and a corresponding movement strategy for the snail. A known result for this
    type of problem states that for an odd number of minutes L, the maximal
    distance is L + (L-1)/2.
    """
    L = 7  # Total duration in minutes

    # The formula for the maximal distance for an odd duration L
    extra_distance = (L - 1) / 2
    max_distance = L + extra_distance

    print("The problem is set over a duration of L minutes.")
    print(f"L = {L}")
    print("\nThe maximal distance can be found using the formula: L + (L-1)/2 for an odd L.")
    print("This formula is derived from finding an optimal arrangement of observers that allows the snail to maximize its travel.")
    print("\nCalculation:")
    print(f"Maximal Distance = {L} + ({L} - 1) / 2")
    print(f"Maximal Distance = {L} + {int(extra_distance)}")
    print(f"Maximal Distance = {int(max_distance)}")

solve_snail_problem()