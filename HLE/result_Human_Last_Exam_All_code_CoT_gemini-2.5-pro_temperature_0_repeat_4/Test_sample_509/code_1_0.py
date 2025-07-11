import math
import numpy as np

def construct_section_for_open_disk(k_points, l):
    """
    Constructs a section for the projection map from conf_l(D) to conf_k(D),
    where D is the open unit disk.

    Args:
        k_points (list of tuples): A list of k points, e.g., [(x1, y1), ...].
        l (int): The desired number of points in the final configuration.

    Returns:
        list of tuples: A list of l points, where the first k are the input points.
    """
    k = len(k_points)
    if l <= k:
        print("Error: l must be greater than k.")
        return k_points

    # Validate that input points are in the open unit disk
    for p in k_points:
        if np.linalg.norm(p) >= 1:
            print(f"Error: Point {p} is not in the open unit disk.")
            return None

    # The construction of a section s: conf_k(D) -> conf_l(D)
    # s(p_1, ..., p_k) = (p_1, ..., p_k, q_1, ..., q_{l-k})
    # We need to define the new points q_j continuously based on p_i.

    # Find the maximum norm among the existing points.
    # This is a continuous function of the input points.
    max_norm = 0.0
    if k > 0:
        max_norm = max(np.linalg.norm(p) for p in k_points)

    # Place the new points on a circle with a radius between max_norm and 1.
    # This ensures the new points are distinct from the old points.
    # The radius is chosen continuously.
    new_radius = (max_norm + 1.0) / 2.0

    l_points = list(k_points)
    num_new_points = l - k

    # Add l-k new points, spaced evenly on the circle of radius new_radius.
    # This ensures the new points are distinct from each other.
    for i in range(num_new_points):
        angle = 2 * math.pi * i / num_new_points
        new_x = new_radius * math.cos(angle)
        new_y = new_radius * math.sin(angle)
        l_points.append((new_x, new_y))

    return l_points

# --- Example Usage ---
# A configuration of k=2 points in the open unit disk
k_config = [(0.1, 0.2), (0.5, -0.3)]
k = len(k_config)
# We want to find a corresponding configuration of l=5 points
l_val = 5

# Apply the section map
l_config = construct_section_for_open_disk(k_config, l_val)

if l_config:
    print(f"Original {k} points: {k_config}")
    print(f"Constructed {l_val} points: {l_config}")

    # Verification
    print("\nVerification:")
    print(f"The first {k} points of the new configuration match the original: {l_config[:k] == k_config}")
    
    newly_added_points = l_config[k:]
    print("Newly added points:")
    for i, p in enumerate(newly_added_points):
        # Using np.isclose for floating point comparison
        print(f"  Point {i+1}: ({p[0]:.4f}, {p[1]:.4f}), Norm: {np.linalg.norm(p):.4f}")
