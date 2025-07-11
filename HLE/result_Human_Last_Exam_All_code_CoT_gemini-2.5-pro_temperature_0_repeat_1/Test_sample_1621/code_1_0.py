def solve_matrix_problem():
    """
    This function solves the problem by finding all natural numbers n
    that satisfy the condition n <= rho(n), where rho(n) is the
    Hurwitz-Radon function.
    """

    def get_rho(n):
        """
        Calculates the value of the Hurwitz-Radon function rho(n).
        n is written as n = 2^k * m, where m is odd.
        k is written as k = 4a + b, where 0 <= b <= 3.
        rho(n) is then 8a + 2^b.
        """
        if n <= 0:
            return 0
        
        # Find k, the exponent of 2 in the prime factorization of n.
        k = 0
        temp_n = n
        while temp_n > 0 and temp_n % 2 == 0:
            k += 1
            temp_n //= 2
            
        # Decompose k into 4a + b.
        a = k // 4
        b = k % 4
        
        # Calculate rho(n).
        return 8 * a + 2**b

    print("The problem is equivalent to finding natural numbers n such that n <= rho(n),")
    print("where rho(n) is the Hurwitz-Radon function.")
    print("We will now test values of n to find the solutions.")
    print("-" * 50)

    solutions = []
    # Based on the mathematical analysis, solutions are expected to be small.
    # We search up to n=16 to demonstrate that the inequality fails for n > 8.
    limit = 16
    
    for n in range(1, limit + 1):
        rho_n = get_rho(n)
        is_solution = n <= rho_n
        
        # The final equation is the inequality n <= rho(n)
        print(f"For n = {n}:")
        print(f"  rho({n}) = {rho_n}")
        print(f"  Checking the inequality: {n} <= {rho_n}")
        
        if is_solution:
            solutions.append(n)
            print(f"  Result: The condition is MET.")
        else:
            print(f"  Result: The condition is NOT MET.")
        print("-" * 20)

    print("\nSummary:")
    print(f"The set of natural numbers n that satisfy the condition is: {solutions}")
    print(f"The total number of such natural numbers is: {len(solutions)}")

solve_matrix_problem()
<<<4>>>