import sys

def solve_knots():
  """
  This function identifies which of the knots K_1 to K_6 are isotopic to the unknot.
  
  The analysis is as follows:
  - K1: This is the 6_1 knot (Stevedore knot), which is not the unknot.
  - K2: This knot can be simplified using Reidemeister moves (type I and type II) to a simple loop. It is the unknot.
  - K3: This knot can be simplified using a Reidemeister type II move. It is the unknot.
  - K4: This is the 4_1 knot (figure-eight knot), which is not the unknot.
  - K5: This knot can be simplified using two Reidemeister type I moves. It is the unknot.
  - K6: This is the 6_2 knot, which is not the unknot.

  Therefore, the indices of the knots isotopic to the unknot are 2, 3, and 5.
  """
  
  # The list of indices of knots that are isotopic to the unknot.
  unknot_indices = [2, 3, 5]
  
  # Print the final list as the answer.
  print(unknot_indices)

solve_knots()