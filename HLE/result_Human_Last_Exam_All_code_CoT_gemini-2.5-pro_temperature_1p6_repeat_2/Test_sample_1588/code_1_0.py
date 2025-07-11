import math

def prime_factor_counts(num):
    """Returns a dictionary of prime factors and their counts for a given number."""
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
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    
    counts = prime_factor_counts(n)
    for p in counts:
        if counts[p] > 1:
            return 0  # n has a squared prime factor
    
    # n is a product of k distinct primes
    k = len(counts)
    return (-1)**k

def get_divisors(n):
    """Returns a sorted list of all divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_problem():
    """
    Calculates the number of non-zero BCH coefficients for a given order
    and prints the step-by-step derivation.
    """
    # The order of the coefficients in the BCH expansion
    n = 10
    # The number of generators (X and Y)
    k = 2

    print(f"To find the number of non-zero coefficients of order {n} in the Baker-Campbell-Hausdorff expansion,")
    print("we use Witt's formula for the number of basic Lie commutators of degree n.")
    print("\nThe formula is: d(n) = (1/n) * sum_{d|n} [mu(d) * k^(n/d)]")
    print(f"For this problem, n = {n} (the order), k = {k} (for operators X and Y), and mu is the Mobius function.")

    divisors = get_divisors(n)
    print(f"\nThe divisors 'd' of n={n} are: {divisors}")
    print("\nLet's calculate each term in the sum:")

    total_sum = 0
    value_parts = []

    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term_value = mu_d * (k**power)
        total_sum += term_value
        
        # Build the string representation for the calculation step
        term_str = f"mu({d})*({k}^({n}/{d}))"
        value_str = f"{mu_d} * {k}^{power}"
        
        print(f"For d={d}: {term_str} = {value_str} = {term_value}")
        
        # Store parts for the final equation printout
        value_parts.append(str(term_value))
        
    result = total_sum // n

    print("\nNow, we plug these values back into Witt's formula:")
    # Building the final equation string
    sum_calc_str = f"d({n}) = (1/{n}) * [ { ' + '.join(value_parts) } ]".replace('+ -', '- ')
    
    print(sum_calc_str)
    print(f"d({n}) = (1/{n}) * [ {total_sum} ]")
    print(f"d({n}) = {result}")

    print(f"\nThus, the number of nonzero coefficients of order {n} is {result}.")

if __name__ == '__main__':
    solve_bch_problem()
