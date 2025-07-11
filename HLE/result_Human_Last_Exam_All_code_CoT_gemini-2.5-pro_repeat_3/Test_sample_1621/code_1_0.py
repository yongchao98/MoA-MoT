import math

def get_b(n):
    """
    For a given n, finds the exponent b in the prime factorization n = m * 2^b, where m is odd.
    """
    if n == 0:
        return 0
    b = 0
    while n > 0 and n % 2 == 0:
        n //= 2
        b += 1
    return b

def rho(n):
    """
    Calculates the Radon-Hurwitz number rho(n).
    n = (2k+1) * 2^b
    b = 4d + c, where c is in {0, 1, 2, 3}
    rho(n) = 2^c + 8d
    """
    if not isinstance(n, int) or n <= 0:
        return 0
    
    b = get_b(n)
    d = b // 4
    c = b % 4
    
    return (2**c) + (8*d)

def solve_problem():
    """
    Finds the number of natural numbers n for which rho(n) >= n.
    """
    
    solution_numbers = []
    # The inequality rho(n) >= n fails for n > 8. We check up to 20 for demonstration.
    limit = 20
    
    print("Finding natural numbers n such that rho(n) >= n:")
    print("-" * 50)
    
    for n in range(1, limit + 1):
        rho_val = rho(n)
        if rho_val >= n:
            solution_numbers.append(n)
            print(f"For n = {n}:")
            print(f"  rho({n}) = {rho_val}")
            print(f"  The condition rho(n) >= n is met ({rho_val} >= {n})")
            print("-" * 50)

    print("\nThe natural numbers n for which the property holds are:", solution_numbers)
    print("\nThe total number of such natural numbers is:", len(solution_numbers))

if __name__ == "__main__":
    solve_problem()