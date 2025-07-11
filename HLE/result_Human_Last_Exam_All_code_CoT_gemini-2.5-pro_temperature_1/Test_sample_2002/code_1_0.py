import math

def is_perfect_square(n):
    """Checks if a non-negative integer n is a perfect square."""
    if n < 0:
        return False
    if n == 0:
        return True
    # Using integer square root for precision and efficiency
    x = math.isqrt(n)
    return x * x == n

def is_fibonacci(n):
    """
    Checks if a number n is a Fibonacci number using Gessel's test.
    A number n is a Fibonacci number if and only if 5*n^2 + 4 or 5*n^2 - 4 is a perfect square.
    """
    if n < 0:
        return False
    # The game starts with T>=2, so we don't need to check for 0 or 1.
    return is_perfect_square(5 * n * n + 4) or is_perfect_square(5 * n * n - 4)

def find_winning_T_for_player2(limit):
    """
    Finds and prints the values of T up to a given limit for which Player 2 has a winning strategy.
    """
    print(f"The problem describes a game known as Fibonacci Nim.")
    print(f"The second player has a winning strategy if and only if the total number of tokens T is a Fibonacci number.")
    print(f"The game requires T >= 2.")
    print(f"\nThe values of T up to {limit} for which the second player wins are:")
    
    winning_values = []
    for T in range(2, limit + 1):
        if is_fibonacci(T):
            winning_values.append(T)
            
    # Print the numbers found
    for val in winning_values:
        print(val)

if __name__ == '__main__':
    # Set a reasonable limit for T to demonstrate the pattern.
    find_winning_T_for_player2(100)