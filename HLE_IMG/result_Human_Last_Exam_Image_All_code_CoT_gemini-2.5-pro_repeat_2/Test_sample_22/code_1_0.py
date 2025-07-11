import sys

def solve_knots():
  """
  This function identifies which of the given knots are isotopic to the unknot.

  Based on visual analysis using Reidemeister moves:
  - K_1: Not the unknot (composite knot).
  - K_2: Is the unknot (can be simplified by undoing twists and unlinked loops).
  - K_3: Is the unknot (a known diagram of the unknot, can be simplified).
  - K_4: Not the unknot (figure-eight knot).
  - K_5: Is the unknot (simplifies to K_3 via a Reidemeister II move).
  - K_6: Not the unknot (trefoil knot).
  
  The indices of the knots that are isotopic to the unknot are therefore 2, 3, and 5.
  """
  
  # The list of indices for knots K_i that are isotopic to the unknot.
  unknot_indices = [2, 3, 5]
  
  # The problem asks to output the numbers in the final equation.
  # The "equation" here is the formation of the result list.
  # We will print the final list as requested by the example format.
  print(unknot_indices)

solve_knots()