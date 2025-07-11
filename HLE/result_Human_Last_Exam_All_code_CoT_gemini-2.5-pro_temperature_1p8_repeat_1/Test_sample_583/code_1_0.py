import math

def f(n):
    """
    Computes f(n) using the formula from Porras-Bejarano and Rodriguez Caballero (2018).
    f(n) is the maximum number of distinct 2-adic valuations of sums of n integers.
    """
    if n <= 0:
        return 0
    if n == 1:
        return 1

    # s2(n) is the number of 1s in the binary representation of n.
    s2_n = bin(n).count('1')

    # k = floor(log2(n))
    k = n.bit_length() - 1

    # Calculate the summation term from the formula:
    # Sum_{i=0 to k-1} (floor(n / 2^(i+1)) - 1) * 2^(i+1)
    # Using integer arithmetic to avoid float precision issues.
    sum_term = 0
    for i in range(k):
        # term = (floor(n / 2**(i+1)) - 1) * 2**(i+1)
        term = ((n >> (i + 1)) - 1) << (i + 1)
        sum_term += term

    # Full formula for f(n) for n >= 2
    return 1 + n - s2_n + sum_term

def main():
    """
    Calculates the limit of f(n) / (n * log2(n)) by testing large values of n.
    """
    print("This script calculates the value of the limit lim_{n->inf} f(n)/(n*log2(n)).")
    print("We will compute the ratio for several large values of n to observe its behavior.")
    
    test_values = [100, 1000, 10000, 100000, 1000000]

    for n in test_values:
        fn_val = f(n)
        log2_n = math.log2(n)
        ratio = fn_val / (n * log2_n)
        print(f"For n = {n}:")
        print(f"  f(n) = {fn_val}")
        print(f"  n * log2(n) = {n * log2_n:.4f}")
        print(f"  Ratio f(n)/(n*log2(n)) = {ratio:.8f}")
        print("-" * 20)

    print("As n increases, the ratio f(n) / (n * log2(n)) approaches 1.")
    # The final answer is the value of the limit.
    final_answer = 1
    print(f"The value of the limit is: {final_answer}")

if __name__ == '__main__':
    main()