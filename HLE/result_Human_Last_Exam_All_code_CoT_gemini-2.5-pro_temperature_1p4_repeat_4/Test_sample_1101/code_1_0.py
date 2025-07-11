import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n)/Gamma(1.5+n-m).
    The Gamma ratio is calculated as a product to avoid large numbers.
    """
    if n < 0:
        return 0
    
    total = 0.0
    for m in range(n + 1):
        # Calculate binomial coefficient C(n, m)
        try:
            binom_coeff = math.comb(n, m)
        except AttributeError: # math.comb is new in Python 3.8
            if m < 0 or m > n:
                binom_coeff = 0
            else:
                binom_coeff = math.factorial(n) // (math.factorial(m) * math.factorial(n-m))

        # Calculate the product term Gamma(1.5+n)/Gamma(1.5+n-m)
        # which is (n+0.5)(n-0.5)...(n+1.5-m)
        product_term = 1.0
        for k in range(m):
            product_term *= (n + 0.5 - k)
            
        term = ((-1)**m) * binom_coeff * product_term
        total += term
        
    return total

def main():
    """
    Calculates the sum for a range of n, and demonstrates the asymptotic behavior.
    """
    print("This script calculates the sum S_n and the ratio |S_n|/n! for n from 0 to 15.")
    print("This ratio should be bounded by a constant if f(n)=n! is the correct complexity function.")
    print("-" * 40)
    print(f"{'n':>3s} {'S_n':>20s} {'|S_n|/n!':>20s}")
    print("-" * 40)
    
    # Asymptotic bound for the ratio |S_n|/n! is e^0.5 / sqrt(pi)
    asymptotic_bound = math.exp(0.5) / math.sqrt(math.pi)

    for n in range(16):
        sn_val = calculate_sum(n)
        try:
            ratio = abs(sn_val) / math.factorial(n)
            print(f"{n:3d} {sn_val:20.6f} {ratio:20.6f}")
        except ValueError: # math.factorial raises for non-integer
            print(f"{n:3d} {sn_val:20.6f} {'-':>20s}")
    
    print("-" * 40)
    print(f"The theoretical maximum for the ratio is approx: {asymptotic_bound:.6f}")
    
    print("\nThe analysis shows that the sum grows asymptotically like n!.")
    print("Therefore, the function f(n) with the lowest complexity that bounds the sum is f(n) = n! or Gamma(n+1).")
    final_function = "n!" # or "Gamma(n+1)"
    
    print("\nThe function f is:")
    print(f"f(n) = {final_function}")

if __name__ == "__main__":
    main()
