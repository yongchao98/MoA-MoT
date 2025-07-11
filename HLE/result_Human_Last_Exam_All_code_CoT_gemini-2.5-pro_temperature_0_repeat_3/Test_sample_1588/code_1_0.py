import math

def get_prime_factorization(num):
    """Returns a dictionary of prime factors and their powers."""
    factors = {}
    d = 2
    temp = num
    while d * d <= temp:
        while temp % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp //= d
        d += 1
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

def mobius(n):
    """Calculates the Mobius function mu(n)."""
    if n == 1:
        return 1
    factors = get_prime_factorization(n)
    for p in factors:
        if factors[p] > 1:
            return 0
    return (-1)**len(factors)

def get_divisors(n):
    """Returns a sorted list of divisors of n."""
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def solve_bch_coefficients():
    """
    Calculates the number of nonzero coefficients of a given order in the
    BCH expansion using Witt's formula and prints the detailed calculation.
    """
    # For this problem, we need the number of coefficients of order 10.
    n = 10  # The order
    k = 2   # The number of generators (X and Y)

    print(f"The number of nonzero coefficients is given by Witt's formula:")
    print(f"L_k(n) = (1/n) * sum_{{d|n}} (mu(d) * k^(n/d))")
    print(f"For n={n} and k={k}:\n")

    divisors = get_divisors(n)
    total_sum = 0
    
    # To store parts of the equation for printing
    line1_parts = []
    line2_parts = []
    line3_parts = []

    for d in divisors:
        mu_d = mobius(d)
        power = n // d
        k_power_val = k**power
        term = mu_d * k_power_val
        total_sum += term
        
        # Build strings for each step of the equation
        line1_parts.append(f"mu({d})*{k}^({n}/{d})")
        line2_parts.append(f"({mu_d})*{k_power_val}")
        line3_parts.append(str(term))

    # Format and print the calculation steps
    line1 = f"L_{k}({n}) = (1/{n}) * [ " + " + ".join(line1_parts) + " ]"
    
    # Replace "+ -" with "-" for cleaner output
    line2_str = " + ".join(line2_parts).replace("+ -", "- ")
    line2 = f"L_{k}({n}) = (1/{n}) * [ {line2_str} ]"
    
    line3_str = " + ".join(line3_parts).replace("+ -", "- ")
    line3 = f"L_{k}({n}) = (1/{n}) * [ {line3_str} ]"
    
    line4 = f"L_{k}({n}) = (1/{n}) * [ {total_sum} ]"
    
    result = total_sum // n
    line5 = f"L_{k}({n}) = {result}"

    print("Calculation:")
    print(line1)
    print(line2)
    print(line3)
    print(line4)
    print(line5)
    
    print(f"\nThe number of nonzero coefficients of order 10 is {result}.")

if __name__ == '__main__':
    solve_bch_coefficients()
<<<99>>>