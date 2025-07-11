import math

def is_perfect_square(n):
    """Checks if a number is a perfect square.
    
    Returns a tuple (is_square, square_root).
    """
    if n < 0:
        return False, 0
    # Add a small epsilon for floating point inaccuracies
    x = int(math.sqrt(n) + 0.5)
    return x * x == n, x

def find_winning_T_for_player2(limit):
    """
    Finds the initial number of tokens T for which the second player has a winning strategy.
    
    According to the theory of Fibonacci Nim, the second player wins if and only if
    the initial number of tokens T is a Fibonacci number.
    The game requires T >= 2 to start.
    
    This function prints the winning values of T up to a given limit and shows the
    mathematical verification for each.
    """
    print(f"The values of T up to {limit} for which the second player wins are the Fibonacci numbers (>= 2).")
    print("A number n is a Fibonacci number if and only if 5*n^2 + 4 or 5*n^2 - 4 is a perfect square.")
    print("---")

    # The game can only start if T >= 2, as the first player must take < T tokens.
    for T in range(2, limit + 1):
        # Check if T is a Fibonacci number using the property 5*T^2 ± 4
        val1 = 5 * T * T + 4
        is_sq_1, root_1 = is_perfect_square(val1)
        
        val2 = 5 * T * T - 4
        is_sq_2, root_2 = is_perfect_square(val2)

        if is_sq_1:
            # Output the number and the equation that proves it's a winning T for Player 2
            print(f"T = {T}: Second player wins. Check: 5 * {T}^2 + 4 = {val1} = {root_1}^2")
        elif is_sq_2:
            # Output the number and the equation that proves it's a winning T for Player 2
            print(f"T = {T}: Second player wins. Check: 5 * {T}^2 - 4 = {val2} = {root_2}^2")

if __name__ == '__main__':
    # We set a limit to demonstrate the principle.
    # The rule applies for all Fibonacci numbers.
    upper_limit = 100
    find_winning_T_for_player2(upper_limit)