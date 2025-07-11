import matplotlib.pyplot as plt
import matplotlib.patches as patches

def draw_sierpinski_components(ax, level, rect_coords):
    """
    Recursively draws the components for the Sierpinski carpet construction.
    - Set A (Sierpinski carpet) is represented by the remaining squares.
    - Set B (complement) is represented by the removed squares.
    - The intersection (A intersect B) contains the boundaries of removed squares.
    """
    if level == 0:
        # At the base level, this square is part of the carpet (Set A)
        ax.add_patch(patches.Rectangle((rect_coords[0], rect_coords[1]), rect_coords[2], rect_coords[3], facecolor='navy', linewidth=0))
        return

    x, y, w, h = rect_coords
    w_third, h_third = w / 3.0, h / 3.0

    for i in range(3):
        for j in range(3):
            sub_x = x + i * w_third
            sub_y = y + j * h_third
            
            if i == 1 and j == 1:
                # This is the central square to be removed (part of Set B's interior)
                removed_square = patches.Rectangle((sub_x, sub_y), w_third, h_third, facecolor='lightgray', linewidth=0)
                ax.add_patch(removed_square)
                # The boundary of the removed square is in the intersection
                intersection_boundary = patches.Rectangle((sub_x, sub_y), w_third, h_third, fill=False, edgecolor='red', linewidth=1.5)
                ax.add_patch(intersection_boundary)
            else:
                # Recurse for the 8 other squares
                draw_sierpinski_components(ax, level - 1, (sub_x, sub_y, w_third, h_third))

def main():
    """
    Main function to set up the plot and demonstrate the construction.
    """
    # The user can change the level to see more components.
    # Level 1 gives 1 component.
    # Level 2 gives 1 + 8 = 9 components.
    # Level 3 gives 1 + 8 + 64 = 73 components.
    # Level N gives (8^N - 1)/7 components.
    level = 3 

    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Set up the plot
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'Sierpinski Carpet Construction: Level {level}')
    
    # Create proxy artists for the legend
    legend_elements = [
        patches.Patch(facecolor='navy', edgecolor='navy', label='Set A (Sierpinski Carpet)'),
        patches.Patch(facecolor='lightgray', edgecolor='lightgray', label='Interior of Set B (Removed Squares)'),
        patches.Patch(fill=False, edgecolor='red', linewidth=1.5, label='Intersection A ∩ B (Boundaries)')
    ]
    ax.legend(handles=legend_elements, loc='upper right')

    # Initial square is the unit square
    initial_rect = (0, 0, 1, 1)
    
    # Start the recursive drawing
    draw_sierpinski_components(ax, level, initial_rect)

    # Count the number of components in the intersection for this level
    num_components = (8**level - 1) // 7 if level > 0 else 0
    
    print("This visualization shows a construction where two closed, connected sets, A and B, whose union is the square, have an intersection with many components.")
    print("Set A is the Sierpinski Carpet (navy).")
    print("Set B is the union of all removed squares (light gray).")
    print("The intersection A ∩ B contains the boundaries of all removed squares (outlined in red).")
    print(f"\nFor level {level}, the number of disconnected components in the intersection is {num_components}.")
    print("As the level increases, this number grows exponentially, demonstrating that there is no largest finite number of components.")

    plt.show()

if __name__ == '__main__':
    main()
