import math

def double_factorial(n):
    """Computes the double factorial n!!."""
    if n < 0:
        return 0
    if n == 0 or n == -1:
        return 1
    res = 1
    for i in range(n, 0, -2):
        res *= i
    return res

def calculate_sum(n):
    """Calculates the sum S(n)."""
    sum_val = 0.0
    
    # Pre-calculate binomial coefficients
    binom_coeffs = [1] * (n + 1)
    for i in range(1, n + 1):
        binom_coeffs[i] = binom_coeffs[i-1] * (n - i + 1) // i
        
    for m in range(n + 1):
        # Calculate Gamma ratio (n+1/2)_m
        gamma_ratio_val = 1.0
        for i in range(m):
            gamma_ratio_val *= (n + 0.5 - i)
            
        term = binom_coeffs[m] * gamma_ratio_val
        
        if m % 2 == 1:
            sum_val -= term
        else:
            sum_val += term
            
    return sum_val

def main():
    """Main function to perform the analysis."""
    print("Analysis of the sum S(n)")
    print("-" * 50)
    print(f"{'n':<5}{'S(n)':<25}{'(2n-1)!!':<20}{'|S(n)|/(2n-1)!!':<20}")
    print("-" * 50)
    
    for n in range(11):
        s_n = calculate_sum(n)
        df_n = double_factorial(2 * n - 1)
        
        ratio = 0.0
        if df_n != 0:
            ratio = abs(s_n) / df_n
            
        print(f"{n:<5}{s_n:<25.6f}{df_n:<20}{ratio:<20.6f}")

if __name__ == "__main__":
    main()
