import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors of num."""
    factors = {}
    d = 2
    temp_num = num
    while d * d <= temp_num:
        while (temp_num % d) == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_num //= d
        d += 1
    if temp_num > 1:
        factors[temp_num] = factors.get(temp_num, 0) + 1
    return factors

def mobius(n):
    """Calculates the Mobius function mu(n)."""
    if n == 0:
        return 0
    if n == 1:
        return 1
    
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
            
    if len(factors) % 2 == 0:
        return 1
    else:
        return -1

def get_divisors(n):
    """Returns a sorted list of all positive divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_bch_coefficients():
    """
    Calculates the number of nonzero coefficients of a given order in the 
    BCH expansion using Witt's formula and prints the detailed calculation.
    """
    n = 10  # Order of the expansion
    k = 2   # Number of generators (X and Y)

    print(f"Calculating the number of nonzero BCH coefficients for order n={n} and k={k} generators.")
    print("Using Witt's formula: L_n(k) = (1/n) * sum_{d|n} (mu(d) * k^(n/d))\n")

    divisors = get_divisors(n)
    print(f"The positive divisors 'd' of n={n} are: {divisors}\n")
    
    total_sum = 0
    sum_terms = []
    
    print("Calculating each term of the summation mu(d) * k^(n/d):")
    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        term = mu_d * (k ** power)
        total_sum += term
        
        # Storing for the final equation printout
        sum_terms.append(str(term))
        
        print(f"For d={d}: mu({d}) * {k}^({n}/{d}) = {mu_d} * {k}^{power} = {term}")

    result = total_sum // n
    
    print("\nPutting it all together:")
    equation_sum = " + ".join(sum_terms).replace("+ -", "- ")
    print(f"L_{n}({k}) = (1/{n}) * ({equation_sum})")
    print(f"L_{n}({k}) = (1/{n}) * ({total_sum})")
    print(f"L_{n}({k}) = {result}")

    print(f"\nThe number of nonzero coefficients of order {n} is {result}.")

if __name__ == "__main__":
    solve_bch_coefficients()