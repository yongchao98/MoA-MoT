import math

def falling_factorial(x, n):
    """Computes the falling factorial [x]_n."""
    if n == 0:
        return 1.0
    res = 1.0
    for i in range(n):
        res *= (x - i)
    return res

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * [n+1/2]_m
    """
    total_sum = 0
    for m in range(n + 1):
        # Calculate binomial coefficient C(n, m)
        comb = math.comb(n, m)
        
        # Calculate falling factorial [n + 1/2]_m
        ff = falling_factorial(n + 0.5, m)
        
        # Add the term to the sum
        term = ((-1)**m) * comb * ff
        total_sum += term
        
    return total_sum

# Calculate and print S_n for n from 0 to 10
for n in range(11):
    result = calculate_sum(n)
    print(f"S({n}) = {result}")
