def find_valid_n():
    """
    Based on the mathematical theory of self-distributive structures,
    an idempotent, commutative, and left self-distributive magma being
    n-cancellable implies that it is medial if and only if n is odd.

    This code will print the first few such values of n.
    """
    print("The values of n for which the implication holds are all positive odd integers.")
    print("The first 10 such values are:")
    
    count = 0
    n = 1
    while count < 10:
        if n % 2 != 0:
            print(n)
            count += 1
        n += 1

find_valid_n()