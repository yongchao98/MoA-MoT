import math

def calculate_f_n(n):
    """
    Calculates the value of f(n) based on the formula:
    f(n) = 1 + sum_{k=1}^{n-1} floor((n-1)/k)
    This formula is a known result from number theory for the maximum number of
    distinct 2-adic valuations of subset sums.
    """
    if n <= 0:
        return 0
    if n == 1:
        return 1
    
    m = n - 1
    # The sum is equivalent to the sum of the divisor function d(j) for j=1 to m.
    total = sum(m // k for k in range(1, m + 1))
    
    return 1 + total

def main():
    """
    The main function to demonstrate the limit calculation.
    """
    print("This script calculates the limit of f(n) / (n * log_2(n)) as n -> infinity.")
    print("The analytical result for this limit is ln(2).")
    print("-" * 50)
    print("Let's verify this numerically for large values of n.")

    # We will test for n being powers of 10.
    for i in range(2, 7):
        n = 10**i
        fn_val = calculate_f_n(n)
        denominator = n * math.log2(n)
        ratio = fn_val / denominator
        print(f"For n = {n:<8}, the ratio is: {ratio:.8f}")

    ln_2 = math.log(2)
    print("-" * 50)
    print(f"The value of ln(2) is approximately: {ln_2:.8f}")
    print("\nAs n increases, the ratio converges to ln(2).")

    # The problem asks for the final equation to be printed.
    # The limit is ln(2). The numbers involved are 2.
    print("\nThe final equation is lim_{n->inf} f(n) / (n * log_2(n)) = ln(2).")
    print("The numbers from the final equation are:")
    print(2)
    print(2)


if __name__ == "__main__":
    main()
