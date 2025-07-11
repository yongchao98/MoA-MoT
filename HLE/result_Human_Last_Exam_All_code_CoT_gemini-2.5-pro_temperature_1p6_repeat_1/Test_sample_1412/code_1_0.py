import math

def get_prime_factorization(num):
    """Returns the prime factorization of a number as a dictionary."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while (temp % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def phi(n):
    """Computes Euler's totient function."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    result = n
    for p in factors:
        result -= result // p
    return result

def tau(n):
    """Computes the number of divisors."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    num_divisors = 1
    for p in factors:
        num_divisors *= (factors[p] + 1)
    return num_divisors

def get_divisors(n):
    """Returns a list of all divisors of a number."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def count_adjustable_graphs(n):
    """
    Computes the number of non-isomorphic connected 3-regular adjustable graphs
    with 2n vertices. This corresponds to counting the number of non-isomorphic
    connected cubic Z2-covers of the cycle Cn.
    """
    
    # Find k, the exponent of 2 in the prime factorization of n
    if n == 0:
        return 0
    k = 0
    temp_n = n
    while temp_n > 0 and temp_n % 2 == 0:
        k += 1
        temp_n //= 2
    
    # Calculate the sum part of the formula
    sum_phi = 0
    divs = get_divisors(n)
    print("Finding divisors d of n={} such that n/d is odd:".format(n))
    for d in divs:
        if (n // d) % 2 == 1:
            phi_d = phi(d)
            print("  - d = {:<4}, n/d = {:<4} (odd), phi(d) = {}".format(d, n//d, phi_d))
            sum_phi += phi_d
    print("Sum of phi(d) for these divisors: {}\n".format(sum_phi))

    # Calculate the tau part of the formula
    # m = n / 2^k
    m = n // (2**k)
    tau_m = tau(m)
    print("Calculating tau(n / 2^k) = tau({}) = {}".format(m, tau_m))

    # Apply the final formula
    count = (sum_phi + tau_m) // 2
    print("\nFinal formula: (sum_phi + tau_m) / 2 = ({} + {}) / 2".format(sum_phi, tau_m))
    return count

if __name__ == '__main__':
    n_vertices = 2000
    n_cycle_len = n_vertices // 2
    result = count_adjustable_graphs(n_cycle_len)
    print("\nThe number of non-isomorphic graphs is {}".format(result))
    print("\n<<<{}>>>".format(result))
