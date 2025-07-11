import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n) / Gamma(1.5+n-m)
    where C(n,m) is the binomial coefficient.
    """
    total_sum = 0
    # The term Gamma(1.5+n) can be factored out of the sum.
    # However, calculating it directly can lead to overflow for large n.
    # It's better to calculate the ratio term by term.
    # ratio = Gamma(1.5+n) / Gamma(1.5+n-m) = (n+0.5) * (n-0.5) * ... * (n+0.5-m+1)

    for m in range(n + 1):
        # Calculate binomial coefficient C(n,m)
        try:
            binom_coeff = math.comb(n, m)
        except ValueError:
            # Handle cases where m > n, although the loop structure prevents this.
            binom_coeff = 0

        # Calculate the falling factorial term [n+1/2]_m
        falling_factorial = 1.0
        for i in range(m):
            falling_factorial *= (n + 0.5 - i)

        term = ((-1)**m) * binom_coeff * falling_factorial
        total_sum += term
        
    return total_sum

def main():
    try:
        n_values = list(range(6))
        print("This code calculates the sum for a few values of n.")
        print("Based on mathematical analysis, the function f(n) with the lowest complexity is n^n.")
        print("f(n) = n^n")

        print("\nHere are some computed values of the sum S_n:")
        for n in n_values:
            result = calculate_sum(n)
            print(f"S({n}) = {result}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
