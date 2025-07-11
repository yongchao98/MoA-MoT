def solve_equation():
    """
    This function finds integer solutions for the number of vertices N based on the equation N(3-m) = 2.
    It iterates through possible integer values for m and checks if N is a positive integer.
    """
    
    solutions = []
    # We are looking for integer solutions for N and m (m>=1)
    # 3-m must be a divisor of 2. Divisors of 2 are -2, -1, 1, 2.
    divisors = [1, 2, -1, -2]
    
    print("Equation to solve for N: N * (3 - m) = 2\n")
    print("Finding positive integer solutions for N (number of vertices) and m (multiplicity, m>=1).")
    
    for divisor in divisors:
        three_minus_m = divisor
        m = 3 - three_minus_m
        if m < 1:
            continue
        
        N = 2 / three_minus_m
        
        # Check if N is a positive integer
        if N > 0 and N == int(N):
            N = int(N)
            solutions.append({'N': N, 'm': m})
            print(f"Found a solution: m = {m}, N = {N}")

    if not solutions:
        print("No valid solutions found.")
        return

    max_N = 0
    for sol in solutions:
        if sol['N'] > max_N:
            max_N = sol['N']

    print(f"\nThe possible numbers of r-vertices are {[sol['N'] for sol in solutions]}.")
    print(f"The maximum number of vertices labelled 'r' is {max_N}.")

solve_equation()