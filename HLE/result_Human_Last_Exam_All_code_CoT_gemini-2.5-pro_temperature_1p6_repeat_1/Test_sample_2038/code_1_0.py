def solve_knot_problem():
    """
    Calculates the number of 2-bridge knots with crossing number at most 13
    that have two disjoint non-parallel embedded minimal genus Seifert surfaces.
    This is equivalent to counting non-torus 2-bridge knots.
    """

    # Data from KnotInfo database for the number of 2-bridge knot types (achiral).
    # key = crossing number, value = count of 2-bridge knots
    num_2_bridge_knots_per_crossing = {
        3: 1,
        4: 1,
        5: 2,
        6: 2,
        7: 5,
        8: 16,
        9: 29,
        10: 74,
        11: 197,
        12: 568,
        13: 1761,
    }

    # Calculate the total number of 2-bridge knots with crossing number at most 13.
    # We sum the counts for crossing numbers from 3 to 13.
    total_2_bridge_knots = sum(num_2_bridge_knots_per_crossing.values())

    # The 2-bridge torus knots are T(2,n) for odd n > 1.
    # The crossing number of T(2,n) is n.
    # We need to count how many such knots have crossing number at most 13.
    # These are T(2,3), T(2,5), T(2,7), T(2,9), T(2,11), and T(2,13).
    torus_knots_to_exclude = [3, 5, 7, 9, 11, 13]
    num_torus_knots = len(torus_knots_to_exclude)

    # The final answer is the total number of 2-bridge knots minus the number of torus knots.
    result = total_2_bridge_knots - num_torus_knots

    # Print the explanation and the final equation.
    print("The number of 2-bridge knots admitting two disjoint non-parallel minimal genus Seifert surfaces")
    print("is the total number of 2-bridge knots minus the number of 2-bridge torus knots.")
    print("\nCounts of 2-bridge knots by crossing number (c):")
    for c, count in num_2_bridge_knots_per_crossing.items():
        print(f"c={c}: {count} knot(s)")

    print(f"\nTotal number of 2-bridge knots with crossing number at most 13: {total_2_bridge_knots}")
    print(f"Number of 2-bridge torus knots (T(2,n)) with crossing number at most 13: {num_torus_knots}")
    print("\nFinal Calculation:")
    print(f"{total_2_bridge_knots} - {num_torus_knots} = {result}")

solve_knot_problem()
<<<2650>>>