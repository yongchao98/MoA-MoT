import math

def solve_circle_grid_problem(radius):
    """
    Calculates the minimal and maximal numbers of grid cells a circle can cross.
    R is the radius of the circle.
    """
    
    # The number of times the circle crosses grid lines is 8*R for a generic position.
    # This corresponds to the maximal number of cells crossed.
    max_cells = 8 * radius
    
    # The minimal number of cells crossed occurs when the circle's center is
    # placed at (I+0.5, J+0.5) for integers I, J.
    # The formula for this specific case is 8*R - 8.
    min_cells = 8 * radius - 8
    
    print(f"The radius of the circle is R = {radius}.")
    print(f"The maximal number of cells crossed is 8 * R = 8 * {radius} = {max_cells}.")
    print(f"The minimal number of cells crossed is 8 * R - 8 = 8 * {radius} - 8 = {min_cells}.")
    print(f"The final answer is ({min_cells}, {max_cells}).")

# Given radius
R = 500
solve_circle_grid_problem(R)
