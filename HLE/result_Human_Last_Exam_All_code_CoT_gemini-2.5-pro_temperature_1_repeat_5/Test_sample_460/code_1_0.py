def solve():
    """
    Calculates the smallest integer u based on the given parameters m and t.
    """
    # Number of items
    m = 4
    # Threshold for the number of agents assigned to an item in O
    t = 20

    # The smallest integer u is given by the formula (m-1)*t.
    u = (m - 1) * t

    # Print the parameters and the calculation step-by-step.
    print(f"Given parameters:")
    print(f"m (number of items) = {m}")
    print(f"t (agent threshold) = {t}")
    print("\nThe smallest integer u is calculated as follows:")
    print(f"u = (m - 1) * t")
    print(f"u = ({m} - 1) * {t}")
    print(f"u = {m - 1} * {t}")
    print(f"u = {u}")

solve()