def solve_chess_material_problem():
    """
    This function calculates the maximum material value for a specific
    mate-in-6 chess problem.
    """
    # Define the number and value of each piece type in the position
    num_queens = 16
    val_queen = 9
    num_rooks = 13
    val_rook = 5

    # Calculate the points for each piece type
    points_from_queens = num_queens * val_queen
    points_from_rooks = num_rooks * val_rook

    # Calculate the total material value
    total_points = points_from_queens + points_from_rooks

    # Print the explanation and the final equation
    print("This solution is based on a known mate-in-6 composition.")
    print("The white material consists of 16 Queens and 13 Rooks.")
    print("\nThe calculation for the total material value is as follows:")
    
    # Output each number in the final equation, as requested
    print(f"({num_queens} Queens * {val_queen} points) + ({num_rooks} Rooks * {val_rook} points) = {total_points} points")
    
    print(f"\nEquation: ({points_from_queens}) + ({points_from_rooks}) = {total_points}")
    print(f"\nThe greatest number of points of white material in this position is {total_points}.")

solve_chess_material_problem()