import math

def solve():
  """
  This function identifies which of the six knots K_i are isotopic to the unknot.

  Analysis of the knots:
  K_1: Appears to be the 6_2 knot, which is non-trivial. Not the unknot.
  K_2: This is a complex diagram of the unknot. It can be simplified by removing the clasps on the left and right sides, which reduces it to the diagram of K_3. Thus, it is the unknot.
  K_3: This diagram can be simplified to a simple loop by a Reidemeister II move. It is the unknot.
  K_4: This is the figure-eight knot (4_1), a non-trivial knot. Not the unknot.
  K_5: This diagram simplifies to the trefoil knot (3_1) after removing the clasp on the right. The trefoil knot is non-trivial. Not the unknot.
  K_6: This is the Stevedore knot (6_1), a non-trivial knot. Not the unknot.

  Therefore, the indices of the knots isotopic to the unknot are 2 and 3.
  """
  unknot_indices = [2, 3]
  print(unknot_indices)

solve()