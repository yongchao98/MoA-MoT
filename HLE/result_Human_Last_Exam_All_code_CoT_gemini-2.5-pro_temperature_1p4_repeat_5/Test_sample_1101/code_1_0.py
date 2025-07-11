import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * (n choose m) * Gamma(3/2+n)/Gamma(3/2+n-m).
    """
    if n < 0:
        return 0

    total_sum = 0
    # Header for the calculation breakdown
    print(f"For n = {n}, the sum S_{n} is calculated as follows:")
    print("S_{n} = ", end="")
    
    for m in range(n + 1):
        # Calculate n choose m
        try:
            binom_coeff = math.comb(n, m)
        except ValueError:
            # math.comb is efficient, but let's have a fallback for older versions.
            if m < 0 or m > n:
                binom_coeff = 0
            else:
                if m == 0 or m == n:
                    binom_coeff = 1
                elif m > n // 2:
                    m = n - m
                
                res = 1
                for i in range(m):
                    res = res * (n - i) // (i + 1)
                binom_coeff = res

        # Calculate Gamma ratio, which is a falling factorial (n+1/2)_m
        # (n+1/2) * (n-1/2) * ... * (n+1/2 - m + 1)
        falling_factorial = 1.0
        for j in range(m):
            falling_factorial *= (n + 0.5 - j)

        term = ((-1)**m) * binom_coeff * falling_factorial
        total_sum += term
        
        # Print the term in the sum equation
        sign = "-" if m % 2 != 0 else "+"
        # For the first term, print without a leading sign if it's positive
        if m == 0:
            if term < 0:
                print(f"- {abs(term):.4f} ", end="")
            else:
                print(f"{term:.4f} ", end="")
        else:
            print(f"{sign} {abs(term):.4f} ", end="")

    print(f"= {total_sum:.4f}")
    return total_sum

# Example: calculate for n=3
calculate_sum(3)

<<<n**(n/2+1)>>>