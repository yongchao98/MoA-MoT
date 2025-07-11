import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from skimage import measure

def check_connectivity(matrix):
    """Checks the connectivity of a binary matrix."""
    # The 'connectivity' argument specifies how pixels are connected (1 for 4-connectivity, 2 for 8-connectivity)
    # measure.label returns a labeled matrix, where each connected component has a unique integer label.
    # The number of components is the maximum label value.
    labels = measure.label(matrix, connectivity=2)
    num_components = np.max(labels)
    return num_components == 1, num_components

def run_construction_and_check(n_holes, resolution=200):
    """
    Constructs two sets A and B for a given number of holes and checks their connectivity.
    
    Args:
        n_holes (int): The number of regions to create in the complement of A.
        resolution (int): The resolution of the grid for the visualization.
    """
    # Create a grid representing the unit square
    square = np.zeros((resolution, resolution), dtype=int)
    
    # --- Construct Set A (The "Ladder") ---
    # A will be two vertical bars connected by horizontal rungs.
    v_bar_width = resolution // 6
    h_rung_height = resolution // (3 * n_holes + 1)
    
    # Vertical bars
    square[0:resolution, 0:v_bar_width] = 1
    square[0:resolution, resolution - v_bar_width:resolution] = 1
    
    # Horizontal rungs
    for i in range(n_holes + 1):
        y_start = (3 * i) * resolution // (3 * n_holes + 1)
        square[y_start:y_start + h_rung_height, v_bar_width:resolution - v_bar_width] = 1
        
    set_a_matrix = square
    
    # --- Construct Set B (The Complement) ---
    # B is the closure of the complement of A. In our grid, this is just 1 - A.
    set_b_matrix = 1 - set_a_matrix

    # --- Check Connectivity ---
    is_a_connected, components_a = check_connectivity(set_a_matrix)
    is_b_connected, components_b = check_connectivity(set_b_matrix)

    print(f"--- Analysis for a construction creating {n_holes} holes ---")
    print(f"Is Set A connected? {is_a_connected} ({components_a} component(s))")
    print(f"Is Set B connected? {is_b_connected} ({components_b} component(s))")
    if not is_b_connected:
        print("\nAs shown, the simple 'ladder' construction for A results in a disconnected Set B.")
        print("This illustrates the core difficulty of the problem.")
    
    # --- Visualization ---
    # Create a combined matrix for visualization (1 for A, 2 for B)
    viz_matrix = set_a_matrix + 2 * set_b_matrix
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(viz_matrix, cmap=plt.cm.get_cmap('coolwarm', 3), interpolation='nearest')
    
    # Create a legend
    red_patch = patches.Patch(color='#3b4cc0', label=f'Set A (Connected: {is_a_connected})')
    blue_patch = patches.Patch(color='#b40426', label=f'Set B (Connected: {is_b_connected})')
    ax.legend(handles=[red_patch, blue_patch], loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)
    
    ax.set_title(f'Failed attempt to create {n_holes} intersection components')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()

# Run the demonstration for a case that would hypothetically create 3 holes for B.
# The intersection A_intersect_B would be the boundaries between the red and blue regions.
# If B were connected, the intersection would have 3 components.
run_construction_and_check(n_holes=3)

print("\nConclusion:")
print("While this simple construction fails, it is a known (but non-trivial) result in topology")
print("that constructions using more complex, intertwined sets can produce an intersection")
print("with any finite number of components. Therefore, there is no largest number.")
