import math

def calculate_sum_term(n, m):
    """
    Calculates the term (n+1/2)(n-1/2)...(n-m+3/2).
    This is equivalent to Gamma(n+1.5) / Gamma(n+1.5-m).
    """
    if m == 0:
        return 1.0
    val = 1.0
    for i in range(m):
        val *= (n + 0.5 - i)
    return val

def calculate_sn(n):
    """
    Calculates the sum S_n.
    """
    total = 0.0
    for m in range(n + 1):
        # Calculate (-1)^m * comb(n,m) * product_term
        term = math.comb(n, m) * calculate_sum_term(n, m)
        if m % 2 == 1:
            total -= term
        else:
            total += term
    return total

def main():
    """
    Analyzes the growth of S_n and finds a bounding function f(n).
    """
    print("n\t S_n\t\t |S_n|/n!\t\t |S_n|/(n+1)!\t |S_n|/(n+2)!")
    print("-" * 70)
    for n in range(11):
        sn = calculate_sn(n)
        
        # Calculate ratios to test bounding functions
        try:
            ratio_n_fact = abs(sn) / math.factorial(n)
        except ZeroDivisionError:
            ratio_n_fact = float('inf')

        try:
            ratio_n_plus_1_fact = abs(sn) / math.factorial(n + 1)
        except ZeroDivisionError:
            ratio_n_plus_1_fact = float('inf')

        try:
            ratio_n_plus_2_fact = abs(sn) / math.factorial(n + 2)
        except ZeroDivisionError:
            ratio_n_plus_2_fact = float('inf')

        print(f"{n}\t {sn:10.4f}\t {ratio_n_fact:10.4f}\t {ratio_n_plus_1_fact:10.4f}\t {ratio_n_plus_2_fact:10.4f}")

if __name__ == "__main__":
    main()