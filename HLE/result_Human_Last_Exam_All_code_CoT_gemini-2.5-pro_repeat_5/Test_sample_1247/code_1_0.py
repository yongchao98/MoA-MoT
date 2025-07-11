def solve():
    """
    This function calculates the number of 1324-avoiding permutations of length n with k inversions,
    for n=333 and k=3.
    """
    n = 333
    
    # The number of 1324-avoiding permutations of length n with 3 inversions
    # for n >= 4 is given by the formula: 2*n^2 - 9*n + 10.
    
    # Coefficients of the polynomial
    coeff_a = 2
    coeff_b = -9
    coeff_c = 10
    
    # Calculate the result
    result = coeff_a * n**2 + coeff_b * n + coeff_c
    
    # Print the equation with the numbers plugged in
    print(f"{coeff_a} * {n}**2 - {abs(coeff_b)} * {n} + {coeff_c} = {result}")

solve()