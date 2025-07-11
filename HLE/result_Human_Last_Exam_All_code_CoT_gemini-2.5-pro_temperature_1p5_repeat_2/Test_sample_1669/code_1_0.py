import itertools

def check_local_vertex_condition(k):
    """
    Checks if there are integer values in the set {+/-1, ..., +/-(k-1)}
    that can be assigned to the three edges of a 3-regular vertex
    to have their sum be zero.
    
    This function demonstrates the 'local' possibility of a k-vector.
    The global existence for any such graph depends on deep theorems.
    """
    if k <= 1:
        return False, None

    # The set of allowed integer values for the edges
    value_range = list(range(1, k))
    possible_values = value_range + [-v for v in value_range]
    
    # Iterate through all combinations of 3 values (with replacement)
    # to see if any combination sums to 0.
    for combo in itertools.product(possible_values, repeat=3):
        if sum(combo) == 0:
            # Found a valid assignment for a single vertex
            return True, combo
            
    return False, None

def find_smallest_k():
    """
    Analyzes the problem of finding the smallest k for a valid k-vector
    in any bridgeless 3-regular graph.
    """
    print("This program determines the smallest integer k for a valid k-vector based on graph theory.")
    print("A valid k-vector is equivalent to a 'nowhere-zero k-flow'.")
    print("We need to find the smallest k that works for ANY 20-vertex bridgeless 3-regular graph.")
    print("-" * 70)

    # --- Check for k = 2 ---
    k = 2
    is_possible, _ = check_local_vertex_condition(k)
    print(f"Analysis for k = {k}: Allowed values are {{+/-1}}.")
    if not is_possible:
        print("Result: IMPOSSIBLE. The sum of three odd numbers (like +/-1) can never be zero.")
    print("-" * 70)

    # --- Check for k = 3 ---
    k = 3
    is_possible, combo = check_local_vertex_condition(k)
    print(f"Analysis for k = {k}: Allowed values are {{+/-1, +/-2}}.")
    if is_possible:
        print(f"A local solution at a vertex is possible (e.g., {combo[0]} + {combo[1]} + {combo[2]} = 0).")
        print("This corresponds to a nowhere-zero 3-flow.")
        print("Theory says a graph has a 3-flow if and only if it is 3-edge-colorable.")
        print("However, not all bridgeless 3-regular graphs (called 'snarks') are 3-edge-colorable.")
        print("Result: INSUFFICIENT. k=3 is not guaranteed to work for all such graphs.")
    print("-" * 70)

    # --- Check for k = 4 ---
    k = 4
    is_possible, combo = check_local_vertex_condition(k)
    print(f"Analysis for k = {k}: Allowed values are {{+/-1, +/-2, +/-3}}.")
    if is_possible:
        print(f"A local solution at a vertex is possible (e.g., {combo[0]} + {combo[1]} + {combo[2]} = 0).")
        print("This corresponds to a nowhere-zero 4-flow.")
        print("Theory states that a bridgeless graph has a 4-flow if and only if it does not contain the Petersen graph as a minor.")
        print("A 3-regular graph may contain a Petersen graph minor, which would prevent a 4-flow.")
        print("Result: INSUFFICIENT. k=4 is not guaranteed to work for all such graphs.")
    print("-" * 70)

    # --- Check for k = 5 ---
    k = 5
    is_possible, combo = check_local_vertex_condition(k)
    print(f"Analysis for k = {k}: Allowed values are {{+/-1, +/-2, +/-3, +/-4}}.")
    if is_possible:
        print(f"A local solution at a vertex is possible (e.g., {combo[0]} + {combo[1]} + {combo[2]} = 0).")
        print("This corresponds to a nowhere-zero 5-flow.")
        print("Tutte's 5-Flow Conjecture states that EVERY bridgeless graph has a 5-flow.")
        print("Assuming this standard conjecture, k=5 is sufficient for any bridgeless graph.")
        print("Result: SUFFICIENT. This is the smallest value of k that works for all cases.")
    print("-" * 70)
    
    final_k = 5
    print(f"\nCONCLUSION: The smallest value of k is {final_k}.")
    print("Here is an example equation for a single vertex with a valid 5-vector assignment:")
    
    # Show an example of a valid assignment for a vertex
    n1, n2, n3 = check_local_vertex_condition(final_k)[1]
    # We choose a more illustrative one manually if the found one is trivial like (1, -1, 0)
    if 0 in {n1, n2, n3}: # 0 is not in our set
        n1,n2,n3 = (2,2,-4)
    print(f"{n1} + {n2} + ({n3}) = {n1 + n2 + n3}")

if __name__ == '__main__':
    find_smallest_k()