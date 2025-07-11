import math
from fractions import Fraction

def find_min_N():
    """
    Finds the minimum positive integer N for which there exists a pair (n1, n2)
    with n1 <= N and n2 <= N such that the group GG_{n1,n2}(r) is infinite.
    This is equivalent to finding the minimum N for which 1/n1 + 1/n2 < 1.
    """
    N = 1
    while True:
        print(f"Testing N = {N}...")
        found_pair = None
        # To find the smallest pair first for cleaner output, we can sort the search
        # but iterating is sufficient as we only need one such pair to exist.
        for n1 in range(1, N + 1):
            for n2 in range(1, N + 1):
                # The condition for an infinite group is 1/n1 + 1/n2 < 1.
                # To avoid floating point issues, we check n1 + n2 < n1 * n2.
                if n1 + n2 < n1 * n2:
                    found_pair = (n1, n2)
                    break
            if found_pair:
                break
        
        if found_pair:
            n1, n2 = found_pair
            print(f"\nFor N = {N}, a solution pair (n1, n2) was found: ({n1}, {n2}).")
            print("This means S(N) is not empty.")
            
            # Show the calculation for the final equation as requested.
            f1 = Fraction(1, n1)
            f2 = Fraction(1, n2)
            s = f1 + f2
            
            print("\nVerifying the condition 1/n1 + 1/n2 < 1:")
            print(f"1/{n1} + 1/{n2} = {f1.numerator}/{f1.denominator} + {f2.numerator}/{f2.denominator} = {s.numerator}/{s.denominator}")
            
            s_float = float(s)
            print(f"The sum is approximately {s_float:.4f}, which is indeed less than 1.")
            
            print(f"\nTherefore, the minimum N for which S(N) is non-empty is {N}.")
            return N
        else:
            print(f"No solution found for N = {N}. S({N}) is empty.")
        
        N += 1

if __name__ == '__main__':
    find_min_N()