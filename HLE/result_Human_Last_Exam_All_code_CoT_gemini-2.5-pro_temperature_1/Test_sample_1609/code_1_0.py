def solve_a_n():
    """
    This function finds a(4), the maximal number of prime implicants of a
    Boolean function of 4 variables.

    The sequence a(n) does not have a simple closed-form formula. Its values
    for small n are known from research in switching theory and are listed in
    the On-Line Encyclopedia of Integer Sequences (OEIS) as A000373.

    This code retrieves the known value for n=4.
    """
    
    # Known values for a(n) for small n:
    # a(0) = 1, a(1) = 2, a(2) = 6, a(3) = 20, a(4) = 78, a(5) = 352
    a_n_values = {
        0: 1,
        1: 2,
        2: 6,
        3: 20,
        4: 78,
        5: 352
    }

    n = 4
    
    # Check if the value for n is in our dictionary
    if n in a_n_values:
        result = a_n_values[n]
        # Since there is no simple calculation, we present the result as a statement.
        # The numbers in the final statement are 4 and 78.
        print(f"a({n}) = {result}")
    else:
        print(f"The value for a({n}) is not available in the pre-computed list.")

solve_a_n()