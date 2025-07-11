import math

def get_b_and_a(n):
    """
    Decomposes n into 2^a * b where b is odd.
    Returns (b, a).
    """
    a = 0
    while n % 2 == 0:
        n //= 2
        a += 1
    b = n
    return b, a

def get_c_and_d(a):
    """
    Decomposes a into 4c + d.
    Returns (c, d).
    """
    c = a // 4
    d = a % 4
    return c, d

def rho(n):
    """
    Calculates the Radon-Hurwitz number for n.
    """
    if n <= 0:
        return 0
    _b, a = get_b_and_a(n)
    c, d = get_c_and_d(a)
    return 8 * c + 2**d

def find_possible_n_values(max_n_to_check):
    """
    Finds the number of natural numbers n for which the condition holds.
    """
    solutions = []
    print("Checking for possible values of n...\n")

    for n in range(1, max_n_to_check + 1):
        # Case n=1 is a solution.
        if n == 1:
            solutions.append(n)
            print(f"n = {n}: This is a known base case solution.")
            continue
            
        # For odd n > 1, no solution exists.
        if n % 2 != 0:
            print(f"n = {n}: No solution exists for odd n > 1.")
            continue

        # For even n, check the condition n <= rho(n).
        rho_n = rho(n)
        is_solution = (n <= rho_n)
        
        print(f"n = {n}: rho({n}) = {rho_n}. Checking if {n} <= {rho_n}: {is_solution}")
        if is_solution:
            solutions.append(n)
    
    return solutions

if __name__ == "__main__":
    # We check up to n=32 to show the pattern holds for larger numbers.
    possible_n = find_possible_n_values(32)
    
    print("\n-------------------------------------------------")
    print(f"The possible values for n are: {possible_n}")
    print(f"The number of natural numbers n for which the condition holds is: {len(possible_n)}")
    print("The final list of numbers is printed above. They are the solutions to the problem.")
