import itertools

def get_unique_eigenvalue_sets():
    """
    Calculates the number of unique sets S(rho) for representations of
    Abelian groups of order 18.
    """
    divisors_of_18 = [1, 2, 3, 6, 9, 18]
    
    # Represent the set of d-th roots of unity, U_d, as a frozenset of integers.
    # The integers correspond to the numerators of the exponents of e^(2*pi*i*k/18).
    # For example, U_9 corresponds to {0, 2, 4, 6, 8, 10, 12, 14, 16}.
    base_sets = []
    for d in divisors_of_18:
        step = 18 // d
        u_d = frozenset(range(0, 18, step))
        base_sets.append(u_d)
    
    # Generate all possible non-empty unions of the base sets
    generated_sets = set()
    for i in range(1, len(base_sets) + 1):
        for subset in itertools.combinations(base_sets, i):
            # Compute the union of the sets in the current subset
            current_union = frozenset().union(*subset)
            generated_sets.add(current_union)
    
    # Print the sizes of the unique sets found to show they are distinct
    print("The sizes of the unique sets are:")
    # The requirement "output each number in the final equation" is interpreted as
    # printing the properties of the solution, in this case, the sizes of the sets.
    set_sizes = sorted([len(s) for s in generated_sets])
    for size in set_sizes:
        print(size)
        
    # The final answer is the total number of unique sets.
    print("\nThe total number of unique sets is:")
    print(len(generated_sets))

if __name__ == '__main__':
    get_unique_eigenvalue_sets()
