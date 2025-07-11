import math

def get_k_and_type(n1, n2):
    """
    Calculates the order k and determines the group type.
    """
    if n1 <= 1 or n2 <= 1:
        return None, "Finite (trivial)"

    # Handle the special case where the composite is a translation
    if n1 == 2 and n2 == 2:
        k_inv = 0
        k = float('inf')
    else:
        # s = 1/n1 + 1/n2 = (n1+n2)/(n1*n2)
        num = n1 + n2
        den = n1 * n2
        common_divisor = math.gcd(num, den)
        # k is the denominator of s in reduced form
        k = den // common_divisor
        k_inv = 1 / k
    
    sum_val = 1/n1 + 1/n2 + k_inv
    
    # Compare with 1, allowing for floating point inaccuracies
    if abs(sum_val - 1) < 1e-9:
        group_type = "Euclidean"
    elif sum_val < 1:
        group_type = "Hyperbolic"
    else:
        group_type = "Spherical"
        
    return k, sum_val, group_type

def find_min_n():
    """
    Finds the minimum N for which a hyperbolic group exists for some (n1, n2) <= N.
    """
    N = 1
    while True:
        N += 1
        print(f"Checking N = {N}...")
        found_hyperbolic_pair_for_N = False
        
        # We only need to check pairs containing N, as smaller pairs were checked before.
        for n1 in range(2, N + 1):
            for n2 in range(n1, N + 1): # n1 <= n2 to avoid duplicates
                if max(n1, n2) != N:
                    continue

                k, sum_val, group_type = get_k_and_type(n1, n2)

                print(f"  Testing pair (n₁, n₂) = ({n1}, {n2}):")
                print(f"    The order of the composite rotation is k = {k}.")
                print(f"    The sum 1/n₁ + 1/n₂ + 1/k is {1/n1} + {1/n2} + {1/k if k != float('inf') else 0} = {sum_val:.4f}")
                print(f"    The underlying group G* is {group_type}.")
                
                if group_type == "Hyperbolic":
                    found_hyperbolic_pair_for_N = True
        
        if found_hyperbolic_pair_for_N:
            print(f"\nFor N={N}, we found a pair, ({n1}, {n2}), that generates a hyperbolic group.")
            print("This is the first N for which S(N) is non-empty under our hypothesis.")
            print(f"The minimum value for N is {N}.")
            return N

if __name__ == '__main__':
    final_N = find_min_n()
    print(f"\nFinal Answer: {final_N}")