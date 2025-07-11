import math

def describe_intersection_components(n):
    """
    This function describes the n components of the intersection A_n an B_n for two
    closed connected sets A_n and B_n whose union is the unit square.

    Args:
        n (int): The desired number of intersection components. Must be a positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for the number of components.")
        return

    print(f"For n={n}, we can construct two closed, connected sets A and B")
    print("whose union is the unit square, and whose intersection A_n ∩ B_n")
    print(f"has exactly {n} connected components.")
    print("\nThese components can be defined as the following disjoint rectangles:")
    print("-" * 60)

    # Define the geometry of the intersection components
    # They will be rectangles of the form [x_start, x_end] x [y_start, y_end]
    x_start = 0.45
    x_end = 0.55

    # Calculate height and spacing for n rectangles
    # Total height available is 1.0. Let's use 80% of it for components and spacing.
    total_space = 0.8
    total_height_of_components = total_space / 2.0
    total_height_of_gaps = total_space / 2.0
    
    component_height = total_height_of_components / n
    gap_height = total_height_of_gaps / (n - 1) if n > 1 else 0

    # Start the first component at y = 0.1
    current_y = 0.1

    for i in range(n):
        y_start = current_y
        y_end = y_start + component_height
        
        # We use round() for cleaner output, not for strict mathematical boundaries
        print(f"Component {i + 1}: Rectangle([x_start={x_start}, x_end={x_end}], "
              f"[y_start={round(y_start, 4)}, y_end={round(y_end, 4)}])")
        
        # Advance y for the next component
        current_y = y_end + gap_height
        
    print("-" * 60)
    print("\nSet A is the union of these n rectangles and a spine on the left, e.g., [0, 0.45] x [0, 1].")
    print("Set B is the union of these n rectangles, a spine on the right, and filler material,")
    print("making both sets connected and their union the whole square.")


# Example: demonstrate the case for 5 components.
# You can change this number to see the construction for any n.
number_of_components = 5
describe_intersection_components(number_of_components)
