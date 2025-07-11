import sys

def find_lexicographically_least_tuple():
    """
    This function finds the lexicographically least tuple (a_1, b_1, ..., a_l, b_l)
    based on the problem's conditions.

    The core logic is to solve the equation derived from the Euler characteristics:
    2 * sum((1-a_i)(1-b_i) for i in 1..l) = l - 1
    where l is minimal and each M(a_i, b_i) is not "full" (i.e., a_i != 1 and b_i != 1).

    The minimal value for l is 3, leading to the equation:
    (1-a_1)(1-b_1) + (1-a_2)(1-b_2) + (1-a_3)(1-b_3) = 1
    """
    
    # We search for the smallest pairs (a,b) first to ensure the final tuple is minimal.
    # A search limit of 10 is sufficient as the minimal solution will involve small numbers.
    search_limit = 10
    
    # Generate a list of (a,b) pairs, sorted lexicographically, where a, b != 1.
    l_pairs = []
    for a in range(search_limit):
        if a == 1:
            continue
        for b in range(search_limit):
            if b == 1:
                continue
            l_pairs.append((a, b))

    # Pre-calculate the x value, x = (1-a)(1-b), for each pair.
    x_vals = {p: (1 - p[0]) * (1 - p[1]) for p in l_pairs}

    # Iterate through combinations of 3 pairs. To ensure the resulting tuple is
    # lexicographically sorted, we iterate such that p1 <= p2 <= p3.
    # The first solution found will be the minimal one.
    for i in range(len(l_pairs)):
        p1 = l_pairs[i]
        for j in range(i, len(l_pairs)):
            p2 = l_pairs[j]
            for k in range(j, len(l_pairs)):
                p3 = l_pairs[k]

                # Check if the pairs solve the equation for l=3.
                if x_vals[p1] + x_vals[p2] + x_vals[p3] == 1:
                    # Found the lexicographically smallest set of pairs.
                    # Flatten the pairs into the final tuple.
                    result_tuple = p1 + p2 + p3
                    
                    # Print the answer as a flat tuple string with no spaces.
                    print(str(result_tuple).replace(" ", ""))

                    # Exit after finding the first (and therefore minimal) solution.
                    return

if __name__ == '__main__':
    find_lexicographically_least_tuple()
