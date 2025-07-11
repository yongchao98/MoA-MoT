import math

def calculate_av_333_3():
    """
    Calculates the number of 1324-avoiding permutations of length 333
    with 3 inversions using the formula: C(n, 3) + C(n-1, 4) for n=333.
    """
    n = 333
    
    # Calculate the first term: C(333, 3)
    term1 = math.comb(n, 3)
    
    # Calculate the second term: C(332, 4)
    term2 = math.comb(n - 1, 4)
    
    # Calculate the sum
    total = term1 + term2
    
    print(f"The number of 1324-avoiding permutations of length {n} with 3 inversions is given by the formula:")
    print(f"av_{n}^3(1324) = C({n}, 3) + C({n-1}, 4)")
    print(f"C({n}, 3) = {term1}")
    print(f"C({n-1}, 4) = {term2}")
    print(f"Total = {term1} + {term2} = {total}")

calculate_av_333_3()