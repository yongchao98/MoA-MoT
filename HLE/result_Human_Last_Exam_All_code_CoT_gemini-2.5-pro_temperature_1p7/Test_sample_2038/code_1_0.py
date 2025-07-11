def count_special_knots():
    """
    This script solves a knot theory problem by finding the number of 2-bridge knots
    with a crossing number up to 13 that admit two disjoint non-parallel minimal
    genus Seifert surfaces.

    The solution relies on the following established theorems in knot theory:
    1. A knot admits two disjoint, non-parallel, minimal genus Seifert surfaces if and
       only if it is either a composite knot or a doubled knot.
    2. All 2-bridge knots are prime, meaning they cannot be composite.
    3. A knot is both a 2-bridge knot and a doubled knot if and only if it is a
       (2,n)-torus knot, where n is an odd integer and |n| > 1.

    Combining these facts, the problem reduces to counting the number of (2,n)-torus knots
    such that their crossing number, which is |n|, is at most 13. Since a knot and its
    mirror image are considered identical, we only need to count for positive n.
    """
    max_crossing_number = 13
    
    # We are looking for odd integers n > 1 such that the crossing number n is <= 13.
    # These n values correspond to the distinct (2,n)-torus knots that satisfy the criteria.
    qualifying_knots = []
    for n in range(3, max_crossing_number + 1, 2):
        qualifying_knots.append(n)

    print("The knots satisfying the given conditions are the (2,n)-torus knots where n is an odd integer, 1 < n <= 13.")
    print("These knots and their crossing numbers (c) are:")
    for n in qualifying_knots:
        print(f"- The (2,{n})-torus knot, with c = {n}")
        
    count = len(qualifying_knots)
    # The final equation requires printing each number that is summed.
    # In this case, we are counting each qualifying knot, so we sum '1' for each.
    equation_components = ["1"] * count
    equation_string = " + ".join(equation_components)
    
    print("\nThe total number of such knots is the sum of the counts for each qualifying crossing number:")
    print(f"Total = {equation_string} = {count}")

count_special_knots()
<<<6>>>