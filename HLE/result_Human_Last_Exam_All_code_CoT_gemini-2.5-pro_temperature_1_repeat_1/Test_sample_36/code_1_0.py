import itertools

def get_divisors(n):
    """Returns all positive divisors of an integer n."""
    divs = set()
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n // i)
    return sorted(list(divs))

def get_canonical_sets(orders):
    """
    Computes the set of all unique eigenvalue sets for a group.
    An eigenvalue set is a union of U_d sets.
    Its canonical form is the frozenset of maximal orders d in the union.
    """
    canonical_forms = set()
    # Generate all non-empty subsets of the orders
    s_list = list(orders)
    subsets = itertools.chain.from_iterable(
        itertools.combinations(s_list, r) for r in range(1, len(s_list) + 1)
    )

    for subset in subsets:
        # Find the maximal elements in the subset under divisibility
        maximals = set()
        for d1 in subset:
            is_maximal = True
            for d2 in subset:
                if d1 != d2 and d2 % d1 == 0:
                    is_maximal = False
                    break
            if is_maximal:
                maximals.add(d1)
        canonical_forms.add(frozenset(maximals))
    return canonical_forms

def main():
    """
    Calculates the number of unique eigenvalue sets for all Abelian groups of order 18.
    """
    # The exponents of the two non-isomorphic Abelian groups of order 18 are 18 and 6.
    exp1 = 18
    exp2 = 6

    # Get the possible orders of character images for each group type
    orders1 = get_divisors(exp1)
    orders2 = get_divisors(exp2)

    # Calculate the set of unique eigenvalue sets for each group type
    sets1 = get_canonical_sets(orders1)
    sets2 = get_canonical_sets(orders2)

    num_sets1 = len(sets1)
    num_sets2 = len(sets2)

    # Calculate the size of the intersection for the final count
    intersection_size = len(sets1.intersection(sets2))

    # The total number of unique sets is |S1 U S2| = |S1| + |S2| - |S1 intersect S2|
    total_unique_sets = num_sets1 + num_sets2 - intersection_size

    print("Calculation Steps:")
    print(f"1. For Abelian groups with exponent 18 (e.g., Z_18), the number of unique eigenvalue sets is {num_sets1}.")
    print(f"2. For Abelian groups with exponent 6 (e.g., Z_2 x Z_3 x Z_3), the number of unique eigenvalue sets is {num_sets2}.")
    print(f"3. The number of sets common to both types is {intersection_size}.")
    print("\nFinal Equation:")
    print(f"Total Unique Sets = (Sets from type 1) + (Sets from type 2) - (Common Sets)")
    print(f"{num_sets1} + {num_sets2} - {intersection_size} = {total_unique_sets}")
    
    print(f"\nThere are {total_unique_sets} unique sets S(rho) intersect D.")
    print(f"<<<{total_unique_sets}>>>")

if __name__ == "__main__":
    main()
