import sys

def find_lexicographically_least_tuple():
    """
    Solves the problem by searching for the lexicographically smallest tuple.

    The problem simplifies to finding three pairs (a_i, b_i) such that:
    1. a_i, b_i are non-negative integers and not equal to 1.
    2. Sum_{i=1 to 3} (1 - a_i)(1 - b_i) = 1.
    3. The combined tuple (a_1, b_1, a_2, b_2, a_3, b_3) is lexicographically minimal.
       To achieve this, we enforce that the pairs are sorted: (a_1, b_1) <= (a_2, b_2) <= (a_3, b_3).
    """

    # Generate a sorted list of possible (pair, C_value)
    # We only need to search small integers for a and b to find the minimal solution.
    search_limit = 5
    possible_pairs = []
    for a in range(search_limit):
        if a == 1:
            continue
        for b in range(search_limit):
            if b == 1:
                continue
            # Enforce a <= b to generate unique pairs like (0,2) but not (2,0),
            # since C(a,b) = C(b,a). The sorted final tuple will be the same.
            if a <= b:
                pair = (a, b)
                c_val = (1 - a) * (1 - b)
                possible_pairs.append({'pair': pair, 'c': c_val})
    
    # Sort primarily by 'a' and secondarily by 'b' to ensure lexicographical order.
    # The generation loop already does this, but an explicit sort is safer.
    possible_pairs.sort(key=lambda x: x['pair'])

    # Search for three pairs that satisfy the C_1 + C_2 + C_3 = 1 equation.
    # The loops are structured to find the lexicographically smallest combination first.
    n = len(possible_pairs)
    for i in range(n):
        for j in range(i, n):
            for k in range(j, n):
                p1 = possible_pairs[i]
                p2 = possible_pairs[j]
                p3 = possible_pairs[k]

                if p1['c'] + p2['c'] + p3['c'] == 1:
                    # Found the first solution, which must be the lexicographically smallest
                    # because of the ordered search.
                    c1, c2, c3 = p1['c'], p2['c'], p3['c']
                    
                    pair1, pair2, pair3 = p1['pair'], p2['pair'], p3['pair']
                    
                    # The pairs are already sorted lexicographically.
                    final_tuple = pair1 + pair2 + pair3
                    
                    print(f"The equation to solve is C1 + C2 + C3 = 1.")
                    print(f"The three pairs (a_i, b_i) are: {pair1}, {pair2}, {pair3}.")
                    print(f"Their corresponding C_i = (1-a_i)(1-b_i) values are: {c1}, {c2}, {c3}.")
                    print(f"Checking the equation: {c1} + {c2} + {c3} = {c1 + c2 + c3}")
                    
                    # Format the final answer as a flat tuple string
                    final_answer_str = str(final_tuple).replace(" ", "")
                    print(f"\nThe lexicographically least tuple is: {final_answer_str}")
                    return

if __name__ == "__main__":
    find_lexicographically_least_tuple()
