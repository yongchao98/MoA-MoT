def radon_hurwitz(n):
    """
    Computes the Radon-Hurwitz number rho(n).
    """
    if n <= 0:
        return 0
    
    # Find m such that n = (odd_part) * 2^m
    odd_part = n
    m = 0
    while odd_part % 2 == 0:
        odd_part //= 2
        m += 1
        
    # Decompose m = 4b + a
    b = m // 4
    a = m % 4
    
    # Calculate rho(n) = 2^a + 8b
    rho = (2**a) + (8 * b)
    return rho

def find_possible_n(limit=100):
    """
    Finds all natural numbers n up to a given limit for which the required matrices exist.
    This is equivalent to finding n where rho(n) >= n.
    """
    solution_n = []
    print("Checking for n from 1 to {}...".format(limit))
    print("-" * 40)
    for n in range(1, limit + 1):
        rho_n = radon_hurwitz(n)
        if rho_n >= n:
            print(f"For n = {n}:")
            print(f"  Radon-Hurwitz number rho({n}) = {rho_n}")
            print(f"  Condition rho(n) >= n is met ({rho_n} >= {n})")
            print(f"  Therefore, n = {n} is a solution.\n")
            solution_n.append(n)

    print("-" * 40)
    print("The possible values for n are:", solution_n)
    print("The number of such natural numbers is:", len(solution_n))

if __name__ == '__main__':
    find_possible_n()
