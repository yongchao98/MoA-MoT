def solve_chess_riddle():
  """
  This function calculates and prints the point value for the minimal
  material required to solve the specified chess problem.

  The solution is a position where a lone black king is stalemated,
  and all squares except the king's square are attacked by white.

  The position used is by J. Mortensen:
  - White: King at g4, Rook at d2, Rook at f2
  - Black: King at h1
  """

  # Define the point value for each piece type involved in the solution
  point_values = {
      'Rook': 5,
  }

  # List the pieces in the minimal solution (excluding the king)
  pieces = ['Rook', 'Rook']

  # Calculate the total points
  total_points = sum(point_values[piece] for piece in pieces)
  
  # Create the equation string
  equation_parts = [str(point_values[piece]) for piece in pieces]
  equation_string = " + ".join(equation_parts)

  # Print the final calculation
  print("The smallest number of points is achieved with two Rooks.")
  print(f"The calculation is: {equation_string} = {total_points}")

solve_chess_riddle()