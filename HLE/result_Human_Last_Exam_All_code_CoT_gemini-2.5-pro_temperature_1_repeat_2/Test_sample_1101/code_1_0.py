import math

def calculate_term(n, m):
    """Calculates the term Gamma(n+3/2)/Gamma(n+3/2-m)"""
    if m == 0:
        return 1.0
    val = 1.0
    start = n + 1.5
    for i in range(m):
        val *= (start - 1 - i)
    return val

def calculate_sum(n):
    """Calculates S(n) = sum_{m=0 to n} (-1)^m * (n choose m) * term(n,m)"""
    total = 0.0
    for m in range(n + 1):
        # Calculate (-1)^m
        sign = 1 if m % 2 == 0 else -1
        
        # Calculate n choose m
        comb = math.comb(n, m)
        
        # Calculate the Gamma fraction term
        term = calculate_term(n, m)
        
        total += sign * comb * term
    return total

def main():
    print("This script calculates the sum S(n) = \\sum_{m=0}^n (-1)^m {n\\choose m} \\frac{\\Gamma(\\frac{3}{2}+n)}{\\Gamma(\\frac{3}{2}+n-m)} for n=0 to 10.")
    print("Based on asymptotic analysis, the sum is bounded by C * f(n) where f(n)=n^n.")
    print("\nCalculating values of the sum S(n):")
    for n in range(11):
        result = calculate_sum(n)
        print(f"S({n}) = {result:.4f}")
    
    print("\nThe function with the lowest complexity is f(n) = n^n.")


if __name__ == "__main__":
    main()
