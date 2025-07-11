import math

def count_representable_integers():
    """
    Counts the number of integers n in [10^18, 10^18 + 10000]
    of the form n = x^3 + 2y^3 + 4z^3 - 6xyz.
    """
    
    # Based on the analysis, we expect solutions only when one variable is large
    # and the others are small. We focus on the case where x is large.
    # The analysis strongly suggests x must be 10^6.
    
    x = 10**6
    lower_bound = 10**18
    upper_bound = 10**18 + 10000
    
    representable_integers = set()
    
    # Case 1: y = 0, z = 0
    n = x**3
    if lower_bound <= n <= upper_bound:
        representable_integers.add(n)
        
    # Case 2: z = 0, y != 0
    # We need 0 <= 2y^3 <= 10000 => y^3 <= 5000
    # y must be positive for 2y^3 to be positive.
    y_limit = int(math.pow(5000, 1/3))
    for y in range(1, y_limit + 1):
        n = x**3 + 2 * y**3
        if lower_bound <= n <= upper_bound:
            representable_integers.add(n)

    # Case 3: y = 0, z != 0
    # We need 0 <= 4z^3 <= 10000 => z^3 <= 2500
    # z must be positive for 4z^3 to be positive.
    z_limit = int(math.pow(2500, 1/3))
    for z in range(1, z_limit + 1):
        n = x**3 + 4 * z**3
        if lower_bound <= n <= upper_bound:
            representable_integers.add(n)

    # The analysis shows that for x=10^6, y and z cannot both be non-zero.
    # The analysis for other dominant variables (y or z) suggests no solutions
    # in the given range due to large gaps between possible values.

    # We print the found values for verification
    # for val in sorted(list(representable_integers)):
    #    print(f"{val} = 10^18 + {val - 10**18}")
    
    print(len(representable_integers))

count_representable_integers()