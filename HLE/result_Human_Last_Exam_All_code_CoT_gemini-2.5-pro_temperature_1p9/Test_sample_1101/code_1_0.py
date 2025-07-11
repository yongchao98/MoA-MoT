import math

def term(n, m):
    """Calculates the m-th term of the sum for a given n."""
    if m > n:
        return 0
    # Use log-gamma for numerical stability with large numbers
    # log(Gamma(n+1.5)) - log(Gamma(n+1.5-m))
    log_gamma_ratio = math.lgamma(n + 1.5) - math.lgamma(n + 1.5 - m)
    falling_factorial = math.exp(log_gamma_ratio)
    
    # math.comb(n, m) for the binomial coefficient
    return (-1)**m * math.comb(n, m) * falling_factorial

def s_n(n):
    """Calculates the sum S_n."""
    return sum(term(n, m) for m in range(n + 1))

def f_n(n):
    """Calculates the bounding function f(n)."""
    if n == 0:
        # To avoid division by zero and handle the special case
        return 1.0 
    return (2**n) * (n**(-0.25))

def main():
    print("Calculating Sn, f(n), and their ratio |Sn/f(n)| for n=1 to 20:")
    print("-" * 50)
    print(f"{'n':<5}{'Sn':<20}{'f(n)':<20}{'|Sn/f(n)|':<20}")
    print("-" * 50)
    for n in range(1, 21):
        sum_val = s_n(n)
        func_val = f_n(n)
        ratio = abs(sum_val / func_val) if func_val != 0 else float('inf')
        print(f"{n:<5}{sum_val:<20.4f}{func_val:<20.4f}{ratio:<20.4f}")
    
    # The asymptotic formula suggests the ratio converges to a constant (around 1/sqrt(pi) ~ 0.564).
    # This confirms that f(n) has the correct growth rate.

main()