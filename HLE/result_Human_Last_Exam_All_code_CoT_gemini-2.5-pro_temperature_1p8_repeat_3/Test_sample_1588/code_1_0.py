import math

def get_divisors(n):
    """
    Finds all positive divisors of an integer n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def mobius(n):
    """
    Calculates the Mobius function mu(n).
    It checks for squared prime factors and counts the number of
    distinct prime factors.
    """
    if n == 1:
        return 1
    
    prime_factors = 0
    i = 2
    num = n
    while i * i <= num:
        if num % i == 0:
            prime_factors += 1
            num //= i
            # Check for squared factor
            if num % i == 0:
                return 0
        i += 1
    
    if num > 1: # The remaining number is a prime factor
        prime_factors += 1
        
    return -1 if prime_factors % 2 != 0 else 1

def solve_bch_coefficient_count():
    """
    Calculates the number of nonzero coefficients for a given order 'n'
    in the BCH expansion using Witt's formula.
    """
    n = 10  # The order of the coefficients
    k = 2   # The number of generators (X, Y)

    divisors = get_divisors(n)
    
    total_sum = 0
    sum_terms_values = []
    
    # Calculate each term in the summation part of Witt's formula
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (k ** power)
        sum_terms_values.append(term)
        total_sum += term
        
    result = total_sum // n
    
    # Format the final equation for printing
    # Example: (1024 - 32 - 4 + 2) / 10 = 99
    equation_sum_part = " + ".join(map(str, sum_terms_values)).replace('+ -', '- ')
    final_equation = f"({equation_sum_part}) / {n} = {result}"
    
    print("The number of nonzero coefficients is calculated using Witt's formula.")
    print("The final calculation is:")
    print(final_equation)
    print(f"\nThus, the number of nonzero coefficients of order 10 is {result}.")

if __name__ == "__main__":
    solve_bch_coefficient_count()
