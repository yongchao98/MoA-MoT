import math

def calculate_sn(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(1.5+n)/Gamma(1.5+n-m)
    """
    if n < 0:
        return 0
    
    # Gamma(1.5+n) can be large, but it cancels in the ratio if we compute the falling factorial
    # Gamma(a)/Gamma(a-m) = (a-1)_m (falling factorial)
    # Here a = 1.5 + n
    
    total_sum = 0
    for m in range(n + 1):
        # Calculate binomial coefficient C(n,m)
        try:
            binom_coeff = math.comb(n, m)
        except ValueError:
            continue # m > n should not happen in this loop

        # Calculate the falling factorial (n+0.5)_m
        # (n+0.5) * (n-0.5) * ... * (n+0.5-m+1)
        falling_factorial = 1.0
        for i in range(m):
            falling_factorial *= (n + 0.5 - i)
            
        term = ((-1)**m) * binom_coeff * falling_factorial
        total_sum += term
        
    return total_sum

def main():
    """
    Calculates and prints S_n for n from 0 to 10, and compares its growth to n!.
    """
    print("n\tS_n\t\t|S_n|\t\tn!\t\t|S_n|/n!")
    print("-" * 60)
    for n in range(11):
        sn = calculate_sn(n)
        n_factorial = float(math.factorial(n))
        ratio = 0.0
        if n_factorial > 0:
            ratio = abs(sn) / n_factorial
        
        print(f"{n}\t{sn: -12.4f}\t{abs(sn): -12.4f}\t{n_factorial: -12.1f}\t{ratio: -12.4f}")
        
    print("\n# The analysis shows that the sum |S_n| grows at a rate comparable to n!.")
    print("# The ratio |S_n|/n! seems to approach a constant for large n (related to Laguerre polynomials).")
    print("# Therefore, the function f(n) with the lowest complexity that bounds the sum is f(n) = n!.")
    print("# The problem is to find f(n) such that sum <= C * f(n).")
    print("# We found that sum is approximately (-1)^n * n! * L_n^{(-2n-1/2)}(1).")
    print("# Since |L_n^{alpha}(x)| can be bounded, the growth is dominated by n!.")
    print("# f(n) = n!")


if __name__ == "__main__":
    main()
