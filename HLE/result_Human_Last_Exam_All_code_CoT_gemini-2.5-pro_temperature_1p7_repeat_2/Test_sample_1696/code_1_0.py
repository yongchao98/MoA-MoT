import math

def calculate_P(n):
    """
    Calculates the value of P(n) using its closed-form formula P(n) = n!.

    Args:
        n (int): A positive integer.

    Returns:
        int: The value of P(n).
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")
    
    # The complex product P(n) simplifies to n!
    result = math.factorial(n)
    
    # Output the formula and the result
    # We create the string for n! = 1 * 2 * ... * n
    factorial_expansion = " * ".join(map(str, range(1, n + 1)))
    if n == 1:
        print(f"P({n}) = {n}! = 1")
    else:
        print(f"P({n}) = {n}! = {factorial_expansion} = {result}")

def main():
    """
    Main function to get user input and call the calculation function.
    """
    try:
        n_str = input("Enter a positive integer n: ")
        n = int(n_str)
        calculate_P(n)
    except ValueError as e:
        print(f"Error: {e}. Please enter a valid positive integer.")

if __name__ == "__main__":
    main()
