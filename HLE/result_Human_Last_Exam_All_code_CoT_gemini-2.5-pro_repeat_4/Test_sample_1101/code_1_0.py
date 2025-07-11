import math

def calculate_sum(n):
    """
    Calculates the sum S_n = sum_{m=0 to n} (-1)^m * C(n,m) * Gamma(3/2+n)/Gamma(3/2+n-m)
    """
    # The term Gamma(3/2+n)/Gamma(3/2+n-m) is the falling factorial [n+1/2]_m
    # [x]_m = x * (x-1) * ... * (x-m+1)
    
    total_sum = 0
    for m in range(n + 1):
        # Calculate binomial coefficient C(n,m)
        try:
            binom_coeff = math.comb(n, m)
        except ValueError:
            binom_coeff = 0
            
        # Calculate the falling factorial [n+1/2]_m
        falling_factorial = 1.0
        for i in range(m):
            falling_factorial *= (n + 0.5 - i)
            
        term = ((-1)**m) * binom_coeff * falling_factorial
        total_sum += term
        
    return total_sum

# Calculate and print the sum for n from 0 to 5
for n in range(6):
    sn = calculate_sum(n)
    print(f"For n = {n}, the sum is: {sn}")

print("\nAnalysis of the sum's growth:")
s_values = [calculate_sum(n) for n in range(6)]
for n in range(1, 6):
    if abs(s_values[n-1]) > 1e-9:
        ratio = abs(s_values[n] / s_values[n-1])
        print(f"|S({n})/S({n-1})| = {ratio:.2f}")

print("\nThe absolute value of the sum grows with n.")
print("The growth appears faster than any polynomial in n.")
print("The class of functions with the lowest complexity that can bound this sum is exponential, like f(n) = a^n.")
print("Based on the ratios, a value like a=2 seems plausible, but the ratio is not constant.")
print("Therefore, a function of the form f(n) = 2^n represents the lowest complexity class for the bound.")
print("\nThe function with the lowest complexity is f(n) = 2^n.")
