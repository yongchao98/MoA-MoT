import math

def get_power_of_2_in_factorization(n):
    """
    Calculates the exponent of 2 in the prime factorization of n.
    For n = odd_part * 2^m, this function returns m.
    """
    if n <= 0:
        raise ValueError("Input must be a positive integer.")
    m = 0
    # While n is even and not zero
    while n % 2 == 0 and n != 0:
        n //= 2
        m += 1
    return m

def radon_hurwitz_number(n):
    """
    Calculates the Radon-Hurwitz number, rho(n).
    If n = (2k+1) * 2^m and m = c + 4d, with 0 <= c <= 3,
    then rho(n) = 2^c + 8d.
    """
    m = get_power_of_2_in_factorization(n)
    c = m % 4
    d = m // 4
    return 2**c + 8 * d

def solve_problem():
    """
    Solves the problem by finding all natural numbers n for which
    the required matrices exist.
    
    This is equivalent to finding n where n <= rho(n), with rho(n)
    being the Radon-Hurwitz number.
    """
    print("Finding natural numbers n for which the condition holds...")
    print("This is equivalent to finding n such that n <= rho(n), where rho(n) is the Radon-Hurwitz number.")
    print("-" * 50)
    
    possible_n = []
    # We only need to check up to a certain limit, as n grows much faster than rho(n).
    # For n > 8, it can be shown that n > rho(n). We check up to 16 to demonstrate this.
    limit = 16
    
    print(f"{'n':>4} | {'rho(n)':>7} | {'n <= rho(n)?':>15}")
    print(f"{'--':>4} | {'-----':>7} | {'------------':>15}")

    for n in range(1, limit + 1):
        rho_n = radon_hurwitz_number(n)
        condition_met = (n <= rho_n)
        if condition_met:
            possible_n.append(n)
        print(f"{n:>4} | {rho_n:>7} | {str(condition_met):>15}")

    print("-" * 50)
    print("The values of n for which such matrices exist are:", possible_n)
    print("The number of such natural numbers n is:", len(possible_n))

# Run the solver
solve_problem()
