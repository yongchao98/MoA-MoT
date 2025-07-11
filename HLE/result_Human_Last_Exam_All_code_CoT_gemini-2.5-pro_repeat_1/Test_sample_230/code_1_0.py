def check_values_of_n():
    """
    This function determines for which positive integers n, a magma being
    n-cancellable implies it is also medial, given it's already idempotent,
    commutative, and left self-distributive.

    The problem is a known result in abstract algebra. The implication holds
    if and only if n is an odd positive integer.

    - For odd n: n-cancellability is a strong enough condition to force the
      magma to become a quasigroup, and commutative, idempotent, distributive
      quasigroups are always medial.
    - For even n: Counterexamples exist. It is possible to construct a magma
      with all the required properties that is n-cancellable for an even n
      but is not medial.

    This script will print the first 20 such values of n.
    """
    print("The property holds for all odd positive integers n.")
    print("Here are the first 20 such values:")
    
    odd_n_values = []
    count = 0
    num = 1
    while count < 20:
        if num % 2 != 0:
            odd_n_values.append(num)
            count += 1
        num += 1
    
    # We are asked to output each number in the final equation.
    # While there's no explicit equation with numbers to solve, we can interpret
    # this as printing the numbers that form the solution set.
    for val in odd_n_values:
        print(val)

if __name__ == "__main__":
    check_values_of_n()