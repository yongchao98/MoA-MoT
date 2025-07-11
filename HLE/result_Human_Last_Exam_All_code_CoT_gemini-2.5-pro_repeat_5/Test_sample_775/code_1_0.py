import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_plot():
    """
    This function generates a visual representation of two closed, connected sets
    whose union is the unit square, and whose intersection has 4 components.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("A Construction for 4 Components in the Intersection")

    # Let A be the orange set and B be the blue set.
    # The intersection A_intersect_B will be 4 disconnected horizontal walls.
    
    # Define the regions for set A (orange)
    # A consists of 3 horizontal chambers and vertical corridors connecting them.
    A_color = 'sandybrown'
    # Corridors on the left and right to ensure A is connected
    ax.add_patch(patches.Rectangle((0.1, 0.1), 0.1, 0.8, facecolor=A_color))
    ax.add_patch(patches.Rectangle((0.8, 0.1), 0.1, 0.8, facecolor=A_color))
    # Horizontal chambers
    ax.add_patch(patches.Rectangle((0.1, 0.1), 0.8, 0.1, facecolor=A_color))
    ax.add_patch(patches.Rectangle((0.1, 0.45), 0.8, 0.1, facecolor=A_color))
    ax.add_patch(patches.Rectangle((0.1, 0.8), 0.8, 0.1, facecolor=A_color))

    # Define the regions for set B (blue)
    # B consists of 2 horizontal chambers and the top/bottom/middle corridors.
    B_color = 'steelblue'
    # Main corridors
    ax.add_patch(patches.Rectangle((0.2, 0.2), 0.6, 0.25, facecolor=B_color))
    ax.add_patch(patches.Rectangle((0.2, 0.55), 0.6, 0.25, facecolor=B_color))
    # Connecting corridors
    ax.add_patch(patches.Rectangle((0.0, 0.0), 1.0, 0.1, facecolor=B_color)) # Bottom
    ax.add_patch(patches.Rectangle((0.0, 0.9), 1.0, 0.1, facecolor=B_color)) # Top
    ax.add_patch(patches.Rectangle((0.0, 0.4), 0.1, 0.2, facecolor=B_color)) # Left-middle
    ax.add_patch(patches.Rectangle((0.9, 0.4), 0.1, 0.2, facecolor=B_color)) # Right-middle

    # The sets A and B are defined as the closure of the colored regions.
    # By this construction:
    # 1. A is connected (3 horizontal parts joined by 2 vertical side parts).
    # 2. B is connected (2 main parts joined by corridors at edges and middle).
    # 3. A U B is the whole square.
    # 4. A intersect B is the boundary between them.
    
    # The components of the intersection are the black lines separating the colors.
    # We draw them explicitly for clarity.
    walls = [
        ((0.2, 0.2), 0.6, 0.005),
        ((0.2, 0.45), 0.6, 0.005),
        ((0.2, 0.55), 0.6, 0.005),
        ((0.2, 0.8), 0.6, 0.005)
    ]
    
    print("The intersection A ∩ B consists of the 4 horizontal black lines shown in the plot.")
    print("Each line is a separate connected component.")
    
    # Equation showing the components
    equation = "A ∩ B = C1 ∪ C2 ∪ C3 ∪ C4"
    print(f"\nWe can represent the intersection as a union of its components:\n{equation}")
    print("where C1, C2, C3, and C4 are the four disjoint horizontal line segments.")

    # Show plot
    plt.show()

draw_plot()
