import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(3/2+n) / Gamma(3/2+n-m)
    """
    if n < 0:
        return 0
    
    # The term Gamma(3/2+n) can be factored out of the sum.
    # The ratio Gamma(a) / Gamma(b) can be calculated as a product.
    # Gamma(a) / Gamma(a-m) = (a-1)*(a-2)*...*(a-m)
    
    total_sum = 0
    for m in range(n + 1):
        # Calculate binomial coefficient C(n, m)
        comb = math.comb(n, m)
        
        # Calculate Gamma(3/2+n) / Gamma(3/2+n-m)
        gamma_ratio = 1.0
        a = 1.5 + n
        for i in range(m):
            gamma_ratio *= (a - 1 - i)
            
        term = ((-1)**m) * comb * gamma_ratio
        total_sum += term
        
    return total_sum

def main():
    """
    Calculates the sum for n=0 to 15 and shows the ratio |S_n|/n!
    """
    print("n\tS_n\t\tn!\t\t|S_n|/n!")
    print("-" * 50)
    for n in range(16):
        sn = calculate_sum(n)
        n_factorial = float(math.factorial(n))
        ratio = abs(sn) / n_factorial if n_factorial > 0 else abs(sn)
        
        print(f"{n}\t{sn:12.4f}\t{n_factorial:12.2e}\t{ratio:8.6f}")

if __name__ == "__main__":
    main()