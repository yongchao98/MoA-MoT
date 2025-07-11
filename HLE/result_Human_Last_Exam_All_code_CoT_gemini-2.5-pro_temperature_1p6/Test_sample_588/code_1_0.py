import itertools

def solve():
    """
    Finds the lexicographically least tuple (a_1, b_1, ..., a_ell, b_ell)
    satisfying the conditions of the problem.
    """
    # We determined the minimal ell is 3.
    # The equation to solve is C1 + C2 + C3 = 1,
    # where C_i = (1 - a_i)(1 - b_i).
    # We search for the lexicographically smallest tuple (a1,b1,a2,b2,a3,b3).

    search_limit = 5 # A small limit for a_i, b_i is sufficient.

    possible_pairs = []
    for a in range(search_limit):
        if a == 1:
            continue
        for b in range(search_limit):
            if b == 1:
                continue
            c = (1 - a) * (1 - b)
            possible_pairs.append({'a': a, 'b': b, 'c': c})

    # The search loops are ordered to find the lexicographically smallest tuple.
    for p1 in possible_pairs:
        for p2 in possible_pairs:
            for p3 in possible_pairs:
                if p1['c'] + p2['c'] + p3['c'] == 1:
                    # We found the first, and therefore smallest, solution.
                    solution_tuple = (p1['a'], p1['b'], p2['a'], p2['b'], p3['a'], p3['b'])
                    
                    print("Found the minimal solution for ell=3.")
                    print(f"The equation is: (1-a1)(1-b1) + (1-a2)(1-b2) + (1-a3)(1-b3) = 1")
                    print("The three pairs are:")
                    print(f"  (a1, b1) = ({p1['a']}, {p1['b']})  --> C1 = (1-{p1['a']})*(1-{p1['b']}) = {p1['c']}")
                    print(f"  (a2, b2) = ({p2['a']}, {p2['b']})  --> C2 = (1-{p2['a']})*(1-{p2['b']}) = {p2['c']}")
                    print(f"  (a3, b3) = ({p3['a']}, {p3['b']})  --> C3 = (1-{p3['a']})*(1-{p3['b']}) = {p3['c']}")
                    print(f"Verification: {p1['c']} + {p2['c']} + {p3['c']} = {p1['c'] + p2['c'] + p3['c']}")
                    print("\nThe lexicographically least tuple is:")
                    # Print the tuple in the required format
                    print(f"({','.join(map(str, solution_tuple))})")
                    return

solve()