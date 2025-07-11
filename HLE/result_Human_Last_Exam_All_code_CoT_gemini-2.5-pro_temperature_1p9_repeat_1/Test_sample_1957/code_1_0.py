import sys

def solve():
    """
    This function solves the problem by finding the minimum value of 1000m + n.
    """
    # The degree of the polynomial det(F) is D = 48.
    # We are given that m * n = D for some positive integers m and n.
    # We want to minimize 1000m + n.
    D = 48

    # Find all divisors of D for possible values of m.
    divisors = []
    for i in range(1, D + 1):
        if D % i == 0:
            divisors.append(i)

    min_val = sys.maxsize
    best_m = 0
    best_n = 0

    # Iterate through all possible values of m
    for m in divisors:
        n = D // m
        val = 1000 * m + n
        if val < min_val:
            min_val = val
            best_m = m
            best_n = n
            
    # The problem asks to output the numbers in the final equation.
    # The optimal choice is m=1, n=48.
    m = best_m
    n = best_n
    result = 1000 * m + n
    
    print(f"The degree of the polynomial is D = {D}.")
    print(f"We want to minimize 1000*m + n, where m*n = {D}.")
    print(f"The optimal choice for m and n is m = {m}, n = {n}.")
    print("The final calculation is:")
    print(f"1000 * {m} + {n} = {result}")

solve()
