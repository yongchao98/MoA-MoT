import math

def is_perfect_square(n):
    """
    Checks if a number is a perfect square.
    Returns a tuple (is_square, root).
    """
    if n < 0:
        return False, 0
    if n == 0:
        return True, 0
    x = int(math.sqrt(n))
    return x * x == n, x

def check_winning_T(T):
    """
    Checks if T is a winning position for the second player and returns the explanation.
    A position T is winning for P2 if and only if T is a Fibonacci number.
    We use the test: T is a Fibonacci number iff 5*T^2 + 4 or 5*T^2 - 4 is a perfect square.
    """
    if not isinstance(T, int) or T <= 0:
        return False, ""

    # Test 5*T^2 + 4
    val1 = 5 * T**2 + 4
    is_sq1, root1 = is_perfect_square(val1)
    if is_sq1:
        explanation = (f"For T = {T}, the second player wins. "
                       f"This is because T is a Fibonacci number, as 5*({T}^2) + 4 = {val1}, which is a perfect square ({root1}^2).")
        return True, explanation

    # Test 5*T^2 - 4
    val2 = 5 * T**2 - 4
    is_sq2, root2 = is_perfect_square(val2)
    if is_sq2:
        explanation = (f"For T = {T}, the second player wins. "
                       f"This is because T is a Fibonacci number, as 5*({T}^2) - 4 = {val2}, which is a perfect square ({root2}^2).")
        return True, explanation
    
    return False, f"For T = {T}, the first player has a winning strategy."

def find_winning_positions(limit):
    """
    Finds and prints all winning T for the second player up to a given limit.
    """
    print(f"The second player has a winning strategy if and only if T is a Fibonacci number.")
    print(f"Testing numbers from 1 to {limit} to find the winning values for T...")
    print("-" * 30)
    
    winning_values = []
    for t_val in range(1, limit + 1):
        is_win, explanation = check_winning_T(t_val)
        if is_win:
            winning_values.append(str(t_val))
            print(explanation)
    
    print("-" * 30)
    print("Summary:")
    print(f"The values of T up to {limit} for which the second player wins are: {', '.join(winning_values)}")

if __name__ == '__main__':
    # Set a limit to check for T
    search_limit = 100
    find_winning_positions(search_limit)