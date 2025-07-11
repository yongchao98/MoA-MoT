import math

def calculate_minimal_edges(d):
    """
    Calculates the minimal number of new edges required to make G' 2-edge-connected,
    given the even integer d.

    Args:
        d (int): An even integer representing a degree, where d >= 2.
    """
    if not isinstance(d, int) or d < 2 or d % 2 != 0:
        print(f"Invalid input for d={d}. 'd' must be an even integer greater than or equal to 2.\n")
        return

    # Step 1: Calculate the maximum number of leaf blocks, p.
    # p = (3d + 2) / 2
    p_numerator = 3 * d + 2
    p_denominator = 2
    p_value = p_numerator / p_denominator

    # Step 2: Calculate the minimal number of edges to add.
    # Edges = ceil(p / 2) = ceil((3d + 2) / 4)
    result_numerator = 3 * d + 2
    result_denominator = 4
    num_edges = math.ceil(result_numerator / result_denominator)

    print(f"For d = {d}:")
    print("The maximum number of leaf blocks (p) in G' is calculated as (3*d + 2) / 2.")
    print(f"  p = (3 * {d} + 2) / {p_denominator} = {p_numerator} / {p_denominator} = {int(p_value)}")
    print("The minimal number of edges to add is ceil(p / 2).")
    print(f"  Edges = ceil({int(p_value)} / 2) = ceil({p_value / 2}) = {num_edges}")
    print("-" * 40)

if __name__ == '__main__':
    print("This script calculates the minimal number of edges for different values of d.")
    print("-" * 40)
    # Demonstrate the calculation for a few valid values of d.
    calculate_minimal_edges(2)
    calculate_minimal_edges(4)
    calculate_minimal_edges(6)
    calculate_minimal_edges(10)
