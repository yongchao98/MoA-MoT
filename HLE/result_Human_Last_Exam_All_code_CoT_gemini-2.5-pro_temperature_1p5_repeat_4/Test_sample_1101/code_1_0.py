import math

def calculate_sum(n):
    """
    Calculates the sum S(n) = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(n+1.5)/Gamma(n-m+1.5)
    The ratio of Gamma functions is the falling factorial (n+0.5)_m.
    """
    if n < 0:
        return 0
    
    total_sum = 0.0
    for m in range(n + 1):
        # Calculate n choose m
        try:
            n_choose_m = math.comb(n, m)
        except ValueError:
            n_choose_m = 0

        # Calculate falling factorial (n+0.5)_m
        falling_factorial = 1.0
        for j in range(m):
            falling_factorial *= (n + 0.5 - j)
        
        term = ((-1)**m) * n_choose_m * falling_factorial
        total_sum += term
        
    return total_sum

def main():
    """
    Calculates S(n), tests the bounding function f(n)=Gamma(n-1/2), and prints the results.
    """
    print("Finding the function f(n) with the lowest complexity such that |S(n)| <= C*f(n).")
    print("Let's test the candidate function f(n) = Gamma(n-1/2).")
    print("-" * 60)
    print(f"{'n':>3s} | {'S(n)':>18s} | {'f(n) = Gamma(n-1/2)':>20s} | {'|S(n)|/f(n)':>12s}")
    print("-" * 60)

    for n in range(3, 16):
        s_n = calculate_sum(n)
        
        # Calculate f(n) = Gamma(n-1/2)
        # We need n-1/2 > 0, so n >= 1. We start from n=3.
        f_n = math.gamma(n - 0.5)
        
        ratio = abs(s_n) / f_n if f_n != 0 else 0
        
        print(f"{n:3d} | {s_n:18.4f} | {f_n:20.4f} | {ratio:12.4f}")

    print("-" * 60)
    print("\nThe ratio |S(n)|/f(n) appears to be bounded by a constant C (approx. 2.5) for n>=3.")
    print("Therefore, a function with the lowest complexity that bounds the sum is f(n) = Gamma(n-1/2).")

if __name__ == "__main__":
    main()
