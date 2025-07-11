import math

def get_prime_factorization(num):
    """
    Returns a dictionary of prime factors of a number.
    e.g., for 12, it returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while temp_num % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
        factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def mobius(n):
    """
    Computes the Mobius function mu(n).
    """
    if n == 1:
        return 1
    
    prime_factors = get_prime_factorization(n)
    
    # If any prime factor has a power greater than 1, mu(n) is 0
    for p in prime_factors:
        if prime_factors[p] > 1:
            return 0
            
    # Otherwise, mu(n) is (-1)^k where k is the number of distinct prime factors
    return (-1) ** len(prime_factors)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of a number.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_order_10():
    """
    Calculates the number of nonzero coefficients of order 10
    in the Baker-Campbell-Hausdorff expansion using Witt's formula.
    """
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    print("The number of nonzero coefficients is calculated using Witt's formula:")
    print(f"L_n(k) = (1/n) * sum(mu(d) * k^(n/d)) for d|n")
    print(f"For this problem, n = {n} and k = {k}.\n")
    
    divisors = get_divisors(n)
    total_sum = 0
    
    equation_parts = []
    value_parts = []
    
    print(f"The divisors of {n} are: {divisors}\n")
    
    print("Calculating terms for each divisor d:")
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_value = mu_d * (k ** power)
        total_sum += term_value
        
        print(f"d = {d}: mu({d}) = {mu_d}, term = {mu_d} * {k}^({n}/{d}) = {mu_d} * {k}^{power} = {term_value}")

        equation_parts.append(f"({mu_d} * {k}^{power})")
        value_parts.append(str(term_value))

    final_result = total_sum // n
    
    print("\nPutting it all together:")
    full_equation = f"L_{n}({k}) = (1/{n}) * [ {' + '.join(equation_parts)} ]"
    print(full_equation.replace('+ -', '- '))
    
    full_values = f"L_{n}({k}) = (1/{n}) * [ {' + '.join(value_parts)} ]"
    print(full_values.replace('+ -', '- '))

    print(f"L_{n}({k}) = (1/{n}) * {total_sum}")
    print(f"L_{n}({k}) = {final_result}")

    print(f"\nThe number of nonzero coefficients of order 10 is: {final_result}")

if __name__ == "__main__":
    solve_bch_order_10()