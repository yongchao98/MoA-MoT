import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n)/Gamma(1.5+n-m)
    """
    if n == 0:
        return 1.0
    
    total_sum = 0.0
    for m in range(n + 1):
        # Binomial coefficient C(n,m)
        try:
            comb = math.comb(n, m)
        except ValueError:
            comb = 0

        # Ratio of Gamma functions as a falling factorial
        # Gamma(A) / Gamma(A-m) = (A-1)...(A-m)
        # Here A = n + 1.5
        prod = 1.0
        for j in range(1, m + 1):
            prod *= (n + 1.5 - j)

        term = ((-1)**m) * comb * prod
        total_sum += term
        
    return total_sum

def f(n):
    """
    Calculates the bounding function f(n) = e^(n/2) * n^(-1/4)
    """
    if n == 0:
        # To avoid division by zero for n^(-1/4)
        # The sum at n=0 is 1. We can define f(0)=1.
        return 1.0
    return math.exp(n / 2.0) * (n**(-0.25))

def main():
    print("Verifying the asymptotic behavior of the sum.")
    print("f(n) = e^(n/2) * n^(-1/4)")
    print("-" * 50)
    print(f"{'n':>3s} | {'Sum S_n':>15s} | {'f(n)':>15s} | {'|S_n| / f(n)':>15s}")
    print("-" * 50)
    
    for n in range(1, 31):
        sum_val = calculate_sum(n)
        f_val = f(n)
        ratio = abs(sum_val) / f_val if f_val != 0 else 0
        
        print(f"{n:3d} | {sum_val:15.4e} | {f_val:15.4e} | {ratio:15.4f}")

if __name__ == "__main__":
    main()
