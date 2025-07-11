import math

def count_cool_strings(n):
    """
    Calculates the number of "cool strings" of maximal length 3n.

    A string is cool if between any two occurrences of the same symbol,
    all symbols are distinct. The number of such strings is n! * 2^n.

    Args:
        n: The number of distinct symbols. Must be a non-negative integer.
    """
    if not isinstance(n, int) or n < 0:
        print("Please enter a non-negative integer for n.")
        return

    # Step 1: Calculate n!
    # There are n! ways to choose the initial permutation of the n symbols.
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        print("n must be a non-negative integer.")
        return

    # Step 2: Calculate 2^n
    # For each of the n symbols, there are 2 choices for its structural
    # arrangement, leading to 2^n possibilities.
    power_of_2 = 2**n

    # Step 3: Multiply the results
    total_strings = n_factorial * power_of_2

    # Output the explanation and the result
    print(f"For n = {n}, the number of cool strings of maximal length is calculated by the formula: n! * 2^n")
    print("The components are:")
    print(f"n! = {n_factorial}")
    print(f"2^n = {power_of_2}")
    print("\nThe final equation is:")
    print(f"{n_factorial} * {power_of_2} = {total_strings}")
    print(f"\nThus, there are {total_strings} cool strings of maximal length for n = {n}.")

if __name__ == '__main__':
    try:
        # Prompt the user for input n
        n_input = int(input("Enter the number of symbols (n): "))
        count_cool_strings(n_input)
    except ValueError:
        print("Invalid input. Please enter an integer.")
