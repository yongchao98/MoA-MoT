import math

def get_radon_hurwitz_number(n):
    """
    Calculates the Radon-Hurwitz number rho(n).
    
    A set of n real n-by-n matrices satisfying the condition exists if and only if
    n is less than or equal to rho(n).

    To calculate rho(n):
    1. Write n = odd_part * 2^b
    2. Write b = 4c + d, where d is in {0, 1, 2, 3}
    3. rho(n) = 2^d + 8c
    """
    if n <= 0:
        return 0
    
    # Step 1: Find b
    b = 0
    temp_n = n
    if temp_n > 0:
        while temp_n % 2 == 0:
            b += 1
            temp_n //= 2
            
    # Step 2: Find c and d from b
    c = b // 4
    d = b % 4
    
    # Step 3: Calculate rho(n)
    rho_n = 2**d + 8 * c
    return rho_n

def solve_problem():
    """
    Finds the number of natural numbers n that satisfy the problem's condition
    by checking the inequality n <= rho(n).
    """
    solutions = []
    print("Checking for which natural numbers n the condition holds...")
    print("The condition is met if and only if n <= rho(n), where rho(n) is the Radon-Hurwitz number.")
    print("-" * 50)

    # The analysis shows that solutions must be small, so a search up to 100 is more than enough
    # to find all of them and be confident in the result.
    for n in range(1, 101):
        rho_n = get_radon_hurwitz_number(n)
        if n <= rho_n:
            solutions.append(n)
            print(f"For n = {n}: rho(n) = {rho_n}. Condition {n} <= {rho_n} is TRUE.")
    
    print("-" * 50)
    print(f"\nThe values of n for which such matrices exist are: {solutions}")
    print(f"In total, there are {len(solutions)} such natural numbers.")

solve_problem()
