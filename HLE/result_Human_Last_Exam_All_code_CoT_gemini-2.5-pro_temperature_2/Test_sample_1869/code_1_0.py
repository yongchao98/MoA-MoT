def find_special_prime():
    """
    Finds the smallest prime p from a list such that Z[p-th root of 6] is not
    the ring of integers of Q(p-th root of 6).

    This is equivalent to finding the smallest prime p in the list that
    satisfies the congruence: 6^(p-1) ≡ 1 (mod p^2).
    """

    choices = [
        ("A", 17),
        ("B", 383),
        ("C", 1093),
        ("D", 66161),
        ("E", 534851)
    ]

    print("The condition for Z[p-th root of 6] NOT being the ring of integers is 6^(p-1) ≡ 1 (mod p^2).")
    print("We will test the given primes in increasing order.\n")

    found_solution = False
    solution_p = -1
    solution_p_squared = -1

    # The choices are already sorted, so the first match will be the smallest prime.
    for label, p in choices:
        p_squared = p * p
        
        # Using Python's built-in pow(base, exp, mod) for modular exponentiation.
        result = pow(6, p - 1, p_squared)
        
        print(f"Testing Choice {label}: p = {p}")
        print(f"Checking if 6^({p-1}) ≡ 1 (mod {p_squared})")
        print(f"Result of 6^{p-1} mod {p_squared} is: {result}")
        
        if result == 1:
            print(f"SUCCESS: p = {p} satisfies the condition.")
            found_solution = True
            solution_p = p
            solution_p_squared = p_squared
            # We break because we need the smallest such prime.
            break 
        else:
            print(f"FAIL: p = {p} does not satisfy the condition.\n")

    if found_solution:
        print("\n----------------------------------------------------")
        print("The smallest prime from the choices that satisfies the condition is", solution_p)
        print("The final equation is:")
        # Final output of the numbers in the equation
        print(f"6^({solution_p - 1}) ≡ 1 (mod {solution_p}^2)")
        print(f"6^{solution_p - 1} ≡ 1 (mod {solution_p_squared})")
        print("----------------------------------------------------")

# Run the function to find the answer
find_special_prime()