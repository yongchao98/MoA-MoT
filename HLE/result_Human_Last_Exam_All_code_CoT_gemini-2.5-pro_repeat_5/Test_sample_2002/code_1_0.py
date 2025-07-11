import math

def is_perfect_square(n):
    """Checks if a number is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    x = int(math.sqrt(n))
    return x * x == n

def is_fibonacci(n):
    """
    Checks if a number is a Fibonacci number using Gessel's test.
    A positive integer n is a Fibonacci number if and only if 
    5*n^2 + 4 or 5*n^2 - 4 is a perfect square.
    The game requires T > 1, so we check for n > 1.
    """
    if n <= 1:
        # T=1 is not a valid start for this game as player 1 must take k < 1 tokens.
        # The first meaningful Fibonacci number for this game is 2.
        return n == 1 
    
    return is_perfect_square(5 * n * n + 4) or is_perfect_square(5 * n * n - 4)

def find_winning_t_for_p2(limit):
    """
    Finds and prints the values of T up to a given limit for which 
    the second player has a winning strategy.
    """
    print(f"Analyzing the game for the total number of tokens T from 2 to {limit}.")
    print("The second player has a winning strategy if and only if T is a Fibonacci number.")
    
    winning_ts = []
    for t in range(2, limit + 1):
        if is_fibonacci(t):
            winning_ts.append(t)
    
    print("\nWithin this range, the values of T for which the second player wins are:")
    # The problem asks to output the numbers, so we print them clearly.
    print(*winning_ts, sep=", ")

if __name__ == '__main__':
    # Set a limit for the total number of tokens T to check.
    T_LIMIT = 150
    find_winning_t_for_p2(T_LIMIT)