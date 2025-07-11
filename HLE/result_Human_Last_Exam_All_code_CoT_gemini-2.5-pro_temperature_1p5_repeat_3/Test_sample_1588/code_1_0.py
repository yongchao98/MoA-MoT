import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num."""
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

def moebius(num):
    """Calculates the Moebius function mu(num)."""
    if not isinstance(num, int) or num < 1:
        raise ValueError("Input must be a positive integer.")
    if num == 1:
        return 1
    
    factors = get_prime_factorization(num)
    for p in factors:
        if factors[p] > 1:
            # Not square-free
            return 0
    
    # Square-free
    if len(factors) % 2 == 1:
        return -1
    else:
        return 1

def get_divisors(num):
    """Returns a sorted list of all positive divisors of num."""
    divs = set()
    for i in range(1, int(math.sqrt(num)) + 1):
        if num % i == 0:
            divs.add(i)
            divs.add(num // i)
    return sorted(list(divs))

def solve_bch_coefficients():
    """
    Calculates the number of nonzero coefficients of a given order in the
    BCH expansion for two generators.
    """
    n = 10  # Order
    k = 2   # Number of generators (X, Y)

    print(f"Calculating the number of nonzero BCH coefficients of order n={n} for k={k} generators.")
    print("Using Witt's formula: L(n, k) = (1/n) * Sum_{d|n} mu(n/d) * k^d\n")

    divisors = get_divisors(n)
    print(f"The positive divisors 'd' of n={n} are: {divisors}")
    print("--------------------------------------------------")
    print("Calculating each term in the summation:")
    
    total_sum = 0
    term_values = []
    for d in divisors:
        mu_val = moebius(n // d)
        k_pow_d = k**d
        term = mu_val * k_pow_d
        total_sum += term
        term_values.append(term)
        
        print(f"For d={d:2}:  mu({n//d:2}) * {k}^{d:<2}  =  {mu_val:2} * {k_pow_d:<4}  =  {term}")

    print("--------------------------------------------------")

    # Format the sum string with proper signs
    sum_str = " + ".join(f"({v})" for v in term_values)
    print(f"The sum of the terms is: {sum_str} = {total_sum}\n")
    
    final_result = total_sum // n
    
    print("Applying the final formula:")
    print(f"L({n}, {k}) = (1/{n}) * {total_sum} = {final_result}")

if __name__ == "__main__":
    solve_bch_coefficients()
    # The final numerical answer required by the system
    # print(f"<<<{99}>>>")