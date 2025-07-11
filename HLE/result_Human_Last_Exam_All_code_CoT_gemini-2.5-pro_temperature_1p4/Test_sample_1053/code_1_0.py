import numpy as np

def calculate_path():
    """
    Calculates and prints the steps to generate a path from a series of movements.
    """
    # Define the series of relative movements (Δx, Δy)
    delta_movements = [
        (0.15, 0.26), (0.38, 0.10), (1.32, -0.06), (0.53, -0.07),
        (1.06, -1.02), (0, 1.32), (-1.06, -1.06), (-0.26, -0.13),
        (-1.12, -0.05), (0.28, -0.48), (0.02, -0.15), (-0.21, 0.26),
        (-0.03, -0.34), (-0.14, 0.26), (-0.17, -0.19), (0.23, 0.38),
        (0.02, 0.25), (-0.78, 0.17), (-0.22, 0.26), (0, 0.26)
    ]

    # Initialize the starting point
    x0, y0 = 0.0, 0.0
    current_x, current_y = x0, y0
    
    print(f"Initial starting point (x_0, y_0) = ({x0:.2f}, {y0:.2f})\n")

    # Sequentially calculate each new point
    for i, (delta_x, delta_y) in enumerate(delta_movements):
        print(f"--- Step {i+1} ---")
        print(f"Current point (x_{i}, y_{i}) = ({current_x:.2f}, {current_y:.2f})")
        print(f"Movement (Δx_{i}, Δy_{i}) = ({delta_x:.2f}, {delta_y:.2f})")
        
        next_x = current_x + delta_x
        next_y = current_y - delta_y
        
        # Print the equations for the next point
        print(f"x_{i+1} = {current_x:.2f} + {delta_x:.2f} = {next_x:.2f}")
        print(f"y_{i+1} = {current_y:.2f} - {delta_y:.2f} = {next_y:.2f}")
        print(f"New point (x_{i+1}, y_{i+1}) = ({next_x:.2f}, {next_y:.2f})\n")
        
        current_x, current_y = next_x, next_y

    # Close the shape by returning to the starting point
    print("--- Closing the Shape ---")
    print(f"Final point (x_20, y_20) = ({current_x:.2f}, {current_y:.2f})")
    print(f"Drawing line segment back to starting point (x_0, y_0) = ({x0:.2f}, {y0:.2f})")

# Execute the function
calculate_path()