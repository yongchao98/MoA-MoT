import math

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    # Find divisors up to the square root of n
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def mobius(n):
    """
    Calculates the Möbius function mu(n).
    - mu(n) = 1 if n is a square-free positive integer with an even number of prime factors.
    - mu(n) = -1 if n is a square-free positive integer with an odd number of prime factors.
    - mu(n) = 0 if n has a squared prime factor.
    """
    if n == 1:
        return 1
    
    prime_factors = []
    d = 2
    temp_n = n
    while d * d <= temp_n:
        if temp_n % d == 0:
            prime_factors.append(d)
            temp_n //= d
            # If the same factor divides again, n is not square-free
            if temp_n % d == 0:
                return 0
        d += 1
    if temp_n > 1:
        prime_factors.append(temp_n)
        
    # Return based on the parity of the number of distinct prime factors
    return -1 if len(prime_factors) % 2 != 0 else 1

def calculate_bch_coefficients(n, k):
    """
    Calculates the number of nonzero coefficients for a given order n and k generators
    using Witt's formula and prints the step-by-step calculation.
    """
    print(f"Calculating the number of nonzero BCH coefficients of order n={n} for k={k} generators.")
    print("Using Witt's formula: L_n(k) = (1/n) * sum over d|n (mu(d) * k^(n/d))\n")

    divisors = get_divisors(n)
    print(f"The divisors 'd' of n={n} are: {divisors}\n")
    
    total_sum = 0
    equation_parts = []
    term_values = []

    print("Calculating each term 'mu(d) * k^(n/d)' for each divisor d:")
    for d in divisors:
        mu_d = mobius(d)
        power_val = n // d
        k_power_n_div_d = k ** power_val
        term = mu_d * k_power_n_div_d
        total_sum += term
        
        print(f"d={d}: mu({d}) * {k}^({n}/{d}) = {mu_d} * {k}^{power_val} = {mu_d} * {k_power_n_div_d} = {term}")
        
        equation_parts.append(f"({mu_d} * {k}^{power_val})")
        term_values.append(str(term))

    result = total_sum // n
    
    # Constructing and printing the final equation
    print("\nFinal Equation:")
    # Show the formula with numbers
    equation_str = f"(1/{n}) * [ {' + '.join(equation_parts)} ]"
    print(equation_str)
    
    # Show the formula with terms calculated
    sum_str = ' + '.join(term_values).replace('+ -', '- ')
    print(f"= (1/{n}) * [ {sum_str} ]")
    
    # Show the result of the sum
    print(f"= (1/{n}) * [ {total_sum} ]")
    
    # Show the final answer
    print(f"= {result}\n")

    print(f"The number of nonzero coefficients of order {n} is {result}.")

# Parameters for the specific problem
order_n = 10
generators_k = 2

calculate_bch_coefficients(order_n, generators_k)