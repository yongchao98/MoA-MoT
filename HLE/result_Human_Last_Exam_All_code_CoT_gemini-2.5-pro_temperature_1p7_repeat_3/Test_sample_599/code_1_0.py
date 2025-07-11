def solve_segmented_number():
    """
    Computes the 50th segmented number.

    The segmented numbers are the positive integers excluding those equal to the
    sum of two or more consecutive smaller terms.
    A number has this property if and only if it is a power of 2.
    The sequence of segmented numbers is 1, 2, 4, 8, ... which is 2^0, 2^1, 2^2, 2^3, ...
    The n-th term in this sequence is 2^(n-1).
    We need to find the 50th term.
    """
    # The term we want to find in the sequence
    n = 50
    
    # The base of the sequence is 2
    base = 2
    
    # The exponent for the n-th term is (n-1)
    exponent = n - 1
    
    # Calculate the 50th segmented number
    result = base ** exponent
    
    # As requested, print the final equation including each number
    print(f"The 50th segmented number is determined by the equation: {base} ** {exponent}")
    print(f"The result is: {result}")

solve_segmented_number()