def solve():
    """
    Calculates the probability that the marble escapes.
    The problem can be modeled as a gambler's ruin problem.
    Let P(n) be the probability of escaping (reaching bin 2025) starting from bin n.
    The absorbing boundaries are at a = 2024 (melting) and b = 2025 (escaping).
    The boundary conditions are:
    P(2024) = 0 (melting means no escape)
    P(2025) = 1 (escaping)

    The recurrence relation for P(n) is given by:
    P(n) = sum_{i != 0} (1/3)**|i| * P(n+i)

    The general solution for a symmetric random walk with zero expected jump distance is a linear function:
    P(n) = A*n + B

    We use the boundary conditions to find the constants A and B.
    1. A * 2025 + B = 1
    2. A * 2024 + B = 0

    Subtracting (2) from (1):
    A * (2025 - 2024) = 1 - 0
    A = 1

    Substituting A=1 into (2):
    1 * 2024 + B = 0
    B = -2024

    So, the solution is P(n) = n - 2024.

    The marble starts at bin 0, so we need to find P(0).
    """

    # Define the positions of the torch and the portal
    torch_bin = 2024
    portal_bin = 2025
    start_bin = 0

    # From the derivation, A = 1 and B = -torch_bin
    A = 1
    B = -torch_bin
    
    # The probability function is P(n) = A*n + B
    # We want to find the probability for the starting bin n=0
    probability = A * start_bin + B

    # Output the equation and the result
    # We write P(0) = 1 * 0 + (-2024)
    print(f"The probability function is P(n) = n - {torch_bin}.")
    print(f"For the starting bin n = {start_bin}, the calculation is:")
    print(f"P({start_bin}) = {A} * {start_bin} + ({B})")
    print(f"Result = {probability}")

solve()
