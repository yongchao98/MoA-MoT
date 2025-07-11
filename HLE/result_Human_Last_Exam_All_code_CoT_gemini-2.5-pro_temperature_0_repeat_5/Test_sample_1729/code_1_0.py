def calculate_pm():
    """
    Calculates the probability Pm for a given positive integer m.

    The problem asks for the probability that a sequence is (i,j)-divisible.
    This condition holds if the 4m indices remaining after removing i and j
    can be partitioned into m arithmetic progressions of length 4.

    Our analysis shows that the number of such pairs (i,j) depends on the
    parity of m:
    - If m is odd, there are 3 favorable pairs.
    - If m is even, there are 5 favorable pairs.

    The total number of ways to choose i and j is C(4m+2, 2).
    """
    try:
        m_str = input("Please enter a positive integer value for m: ")
        m = int(m_str)
        if m <= 0:
            print("Error: m must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Determine the number of favorable pairs (i,j) based on the parity of m.
    if m % 2 == 1:
        # m is odd
        numerator = 3
    else:
        # m is even
        numerator = 5

    # Calculate the total number of pairs (i,j), which is C(4m+2, 2).
    denominator = (2 * m + 1) * (4 * m + 1)

    # Output the results, showing each number in the final equation.
    print(f"\nFor m = {m}:")
    print(f"The number of favorable pairs (i, j) is {numerator}.")
    print(f"The total number of pairs (i, j) is (2 * {m} + 1) * (4 * {m} + 1) = {denominator}.")
    print(f"The probability P_{m} is calculated as:")
    print(f"P_{m} = {numerator} / {denominator}")

# Execute the function
calculate_pm()