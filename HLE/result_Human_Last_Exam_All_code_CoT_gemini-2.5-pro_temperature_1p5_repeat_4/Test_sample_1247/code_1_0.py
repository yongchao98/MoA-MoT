import math

def solve():
    """
    Calculates av_333^3(1324) based on established combinatorial theorems and formulas.
    """
    n = 333
    k = 3

    print(f"We want to find av_{n}^{k}(1324).")
    print("This problem can be solved by applying known results from the field of permutation patterns.")
    print(f"1. For k < n, the number of 1324-avoiding permutations of length n with k inversions, av_n^k(1324), equals the number of 123-avoiding permutations, N(n, k).")
    print(f"2. The number of 123-avoiding permutations, N(n, k), equals the number of 132-avoiding permutations, T(n, k).")
    print(f"3. A formula exists for T(n, k) when k=3: T(n, 3) = C(n-1, 3) + n - 2.")
    print(f"\nCombining these, we need to calculate T({n}, {k}) using the formula:")
    
    n_minus_1 = n - 1
    n_minus_2 = n - 2
    
    print(f"T({n}, {k}) = C({n} - 1, {k}) + {n} - 2 = C({n_minus_1}, {k}) + {n_minus_2}")
    
    # Step 1: Calculate the combination C(n-1, k)
    num1 = n_minus_1
    num2 = n_minus_1 - 1
    num3 = n_minus_1 - 2
    
    den1 = 3
    den2 = 2
    den3 = 1
    
    numerator = num1 * num2 * num3
    denominator = den1 * den2 * den3
    
    comb_val = math.comb(n_minus_1, k)

    print(f"\nFirst, we calculate C({n_minus_1}, {k}):")
    print(f"C({n_minus_1}, {k}) = ({num1} * {num2} * {num3}) / ({den1} * {den2} * {den3})")
    print(f"C({n_minus_1}, {k}) = {numerator} / {denominator}")
    print(f"C({n_minus_1}, {k}) = {comb_val}")

    # Step 2: Add the second part of the formula
    result = comb_val + n_minus_2
    
    print(f"\nNow we complete the calculation:")
    print(f"T({n}, {k}) = {comb_val} + {n_minus_2}")
    print(f"The final result is: {result}")

solve()
<<<6044391>>>