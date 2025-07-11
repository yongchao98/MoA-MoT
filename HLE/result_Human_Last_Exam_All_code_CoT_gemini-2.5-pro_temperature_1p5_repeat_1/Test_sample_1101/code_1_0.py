import math

def falling_factorial(x, n):
    """Computes the falling factorial x_n = x * (x-1) * ... * (x-n+1)"""
    if n == 0:
        return 1.0
    res = 1.0
    for i in range(n):
        res *= (x - i)
    return res

def nCr_times_falling_factorial(n, m):
    """
    Computes binomial(n, m) * (n + 1/2)_m without using large intermediate numbers
    for factorials.
    binom(n,m) = n! / (m! * (n-m)!)
    (n+1/2)_m = (n+1/2)*(n-1/2)*...*(n+1/2-m+1)
    
    The term is = (n/m) * ((n-1)/(m-1)) * ... * ((n-m+1)/1) * (n+1/2)_m
    """
    term = 1.0
    # Calculate binomial coefficient part
    for i in range(m):
        term = term * (n - i) / (i + 1)
    # Calculate falling factorial part
    term *= falling_factorial(n + 0.5, m)
    return term
    
def calculate_sum(n):
    """
    Calculates S_n = sum_{m=0 to n} (-1)^m * C(n, m) * (n+1/2)_m
    """
    total = 0.0
    for m in range(n + 1):
        sign = 1 if m % 2 == 0 else -1
        term = nCr_times_falling_factorial(n, m)
        total += sign * term
    return total

print("n\tS_n\t\t|S_n/S_{n-1}|\t\t|S_n|/n!")
s_prev = 1.0 # S_0
for n in range(11):
    s_curr = calculate_sum(n)
    
    # Calculate ratio of absolute values
    ratio = "N/A"
    if n > 0:
        if abs(s_prev) > 1e-9: # Avoid division by zero
            ratio = f"{abs(s_curr/s_prev):.4f}"
        else:
            ratio = "Inf"
    
    # Calculate ratio with n!
    n_factorial = math.factorial(n)
    ratio_factorial = "N/A"
    if n_factorial > 0:
        ratio_factorial = f"{abs(s_curr / n_factorial):.4f}"

    print(f"{n}\t{s_curr:<15.4f}\t{ratio:<15}\t{ratio_factorial:<15}")
    s_prev = s_curr
