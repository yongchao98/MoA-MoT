def solve():
    """
    Solves the problem by finding how many natural numbers n satisfy the given matrix condition.
    This is equivalent to finding n such that the Hurwitz-Radon function rho(n) >= n.
    """

    def rho(n_val):
        """Calculates the Hurwitz-Radon function rho(n)."""
        if n_val <= 0:
            return 0
        
        # Step 1: Decompose n into odd_part * 2^b
        odd_part = n_val
        b = 0
        while odd_part % 2 == 0:
            odd_part //= 2
            b += 1
            
        # Step 2: Decompose b into 4c + d
        c = b // 4
        d = b % 4
        
        # Step 3: Apply the Hurwitz-Radon formula
        return 8 * c + 2**d

    solutions = []
    # We test values of n up to 32 to demonstrate the pattern.
    # The mathematical proof shows there are no solutions beyond n=8.
    limit = 32
    
    print("This program checks for which natural numbers n the condition rho(n) >= n holds.")
    print("The final equation we are testing is `rho(n) >= n`.")
    print("-" * 50)
    print("n\t | rho(n)\t | Condition (rho(n) >= n)")
    print("-" * 50)

    for n in range(1, limit + 1):
        rho_val = rho(n)
        is_solution = rho_val >= n
        if is_solution:
            solutions.append(n)
        
        print(f"{n:<8d} | {rho_val:<12d} | {is_solution}")

    print("-" * 50)
    print(f"\nThe set of solutions is: {solutions}")
    print(f"There are {len(solutions)} such natural numbers.")

solve()