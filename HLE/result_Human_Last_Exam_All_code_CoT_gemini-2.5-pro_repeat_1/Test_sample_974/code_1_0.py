import sympy

def find_prime_divisors():
    """
    Finds all prime divisors p of q=12740347 such that the number of elements
    of order p in PSL(3, q^2) and PSL(4, q) are equal.
    """
    q = 12740347
    print(f"The given number is q = {q}.")

    # Step 1: Find the prime divisors of q.
    # We interpret the question as asking for prime divisors `p` of `q`.
    prime_divisors_of_q = list(sympy.factorint(q).keys())

    print(f"The prime divisors of q are: {prime_divisors_of_q}")

    # Step 2: For each prime divisor p, check the conditions.
    # If p | q, then p is the characteristic of the field F_q.
    # The number of elements of order p in PSL(n, r) (r is a power of q)
    # is r^(n^2-n) - 1, provided that p >= n.

    solutions = []
    for p in prime_divisors_of_q:
        print(f"\n--- Checking for p = {p} ---")

        # Define n for both groups
        n1 = 3
        n2 = 4

        # Check condition for PSL(3, q^2)
        cond1 = (p >= n1)
        print(f"For PSL(3, q^2), n={n1}. The condition for the formula is p >= {n1}.")
        print(f"Checking if {p} >= {n1}: {cond1}")

        # Check condition for PSL(4, q)
        cond2 = (p >= n2)
        print(f"For PSL(4, q), n={n2}. The condition for the formula is p >= {n2}.")
        print(f"Checking if {p} >= {n2}: {cond2}")

        if cond1 and cond2:
            print("\nBoth conditions are met. We can apply the formulas.")
            
            # The number of elements of order p in PSL(3, q^2)
            expr1_n = n1
            expr1_exp = f"{expr1_n}^2 - {expr1_n}"
            expr1_val = expr1_n**2 - expr1_n
            
            # The number of elements of order p in PSL(4, q)
            expr2_n = n2
            expr2_exp = f"{expr2_n}^2 - {expr2_n}"
            expr2_val = expr2_n**2 - expr2_n

            print(f"Number of elements of order p in PSL({n1}, q^2) is (q^2)^({expr1_exp}) - 1 = q^{2*expr1_val} - 1 = q^12 - 1.")
            print(f"Number of elements of order p in PSL({n2}, q) is q^({expr2_exp}) - 1 = q^{expr2_val} - 1 = q^12 - 1.")
            
            print(f"The numbers of elements are equal, so p = {p} is a solution.")
            solutions.append(p)
        else:
            print(f"\nOne or both conditions are not met, so this method does not apply for p = {p}.")

    print("\n--- Final Result ---")
    if not solutions:
        print("There are no such prime divisors.")
    else:
        print(f"The list of all prime divisors satisfying the condition is: {solutions}")
        
        # As requested, outputting the numbers in the final equation
        p_sol = solutions[0]
        print(f"\nFor the solution p = {p_sol}, the final equality is:")
        print(f"Number_of_elements(PSL(3, {q}^2)) = Number_of_elements(PSL(4, {q}))")
        print(f"({q}^2)^({3}^2 - {3}) - 1 = {q}^({4}^2 - {4}) - 1")
        print(f"({q}^2)^{6} - 1 = {q}^{12} - 1")
        print(f"{q}^{12} - 1 = {q}^{12} - 1")


find_prime_divisors()