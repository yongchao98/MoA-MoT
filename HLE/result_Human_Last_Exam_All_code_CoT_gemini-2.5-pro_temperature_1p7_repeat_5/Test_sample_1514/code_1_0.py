def solve():
    """
    This function calculates the smallest number of topologically distinct
    compactifications of the ray with a remainder X.

    This number is equivalent to the number of orbits of non-empty continua
    in X under the action of the homeomorphism group of X. To find the smallest
    number, we choose the simplest possible space X that satisfies the criteria:
    the two-point discrete space X = {0, 1}.
    """

    # 1. The non-empty continua in X = {0, 1} are the singletons {0} and {1}.
    # We use frozensets because they are hashable and can be stored in a set.
    continua = {frozenset({0}), frozenset({1})}
    num_continua = len(continua)

    # 2. The homeomorphisms of X are the identity and the swap map.
    # We represent these as dictionaries mapping points to points.
    h_id = {0: 0, 1: 1}
    h_swap = {0: 1, 1: 0}
    homeomorphisms = [h_id, h_swap]

    # 3. This helper function applies a homeomorphism to a continuum.
    def apply_homeomorphism(h, c):
        return frozenset({h[p] for p in c})

    # 4. We calculate the number of orbits.
    visited = set()
    num_orbits = 0

    print("Chosen space X = {0, 1}")
    print(f"Number of non-empty continua in X: {num_continua}")

    for c in continua:
        if c not in visited:
            num_orbits += 1
            # Found a new orbit. Now, find all continua in this orbit.
            orbit_of_c = set()
            queue = [c]
            processed_in_queue = set([c])

            while queue:
                current_c = queue.pop(0)
                orbit_of_c.add(current_c)
                for h in homeomorphisms:
                    new_c = apply_homeomorphism(h, current_c)
                    if new_c not in processed_in_queue:
                        queue.append(new_c)
                        processed_in_queue.add(new_c)

            # Mark all continua in the found orbit as visited.
            visited.update(orbit_of_c)

    print(f"The number of orbits of continua is: {num_orbits}")
    print("\nSince the number of orbits determines the number of distinct")
    print("compactifications, and we have shown that 1 is achievable and")
    print("is the minimum possible number, the final answer is 1.")

solve()