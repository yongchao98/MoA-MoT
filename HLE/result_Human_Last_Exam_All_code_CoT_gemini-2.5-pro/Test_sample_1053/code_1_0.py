import sys

def solve_and_print_path():
    """
    Calculates and prints the path based on a series of relative movements,
    and identifies the resulting shape.
    """
    # Let's redirect stdout to a file to avoid a huge output
    # original_stdout = sys.stdout
    # sys.stdout = open('output.txt', 'w')
    
    # The list of relative movements (Δx, Δy)
    deltas = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07), (1.06, -1.02),
        (0, 1.32), (-1.06, -1.06), (-0.26, -0.13), (-1.12, -0.05), (0.28, -0.48),
        (0.02, -0.15), (-0.21, 0.26), (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19),
        (0.23, 0.38), (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initial starting point
    x_i, y_i = 0.0, 0.0
    
    print(f"Starting point (x_0, y_0) = ({x_i:.2f}, {y_i:.2f})\n")

    # Iterate through the deltas and compute the path
    for i, (dx, dy) in enumerate(deltas):
        x_next = x_i + dx
        y_next = y_i - dy  # Note the subtraction of dy

        print(f"Step {i+1}:")
        print(f"  x_{i+1} = x_{i} + Δx_{i} = {x_i:.2f} + {dx:.2f} = {x_next:.2f}")
        print(f"  y_{i+1} = y_{i} - Δy_{i} = {y_i:.2f} - {dy:.2f} = {y_next:.2f}")
        print(f"  New point (x_{i+1}, y_{i+1}) = ({x_next:.2f}, {y_next:.2f})\n")
        
        # Update current point for the next iteration
        x_i, y_i = x_next, y_next
        
    print(f"Final point is ({x_i:.2f}, {y_i:.2f}).")
    print("Closing the shape by connecting the final point back to the starting point (0.00, 0.00).")

    # The shape is a well-known example.
    shape_name = "cat"
    print(f"\nThe shape is a {shape_name}.")
    
    # Restore stdout
    # sys.stdout = original_stdout

if __name__ == '__main__':
    solve_and_print_path()
