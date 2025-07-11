import math

def solve():
  """
  This function determines the largest value of I_G.
  Based on the analysis, the group G = Q/Z (rational numbers modulo 1)
  is a countably infinite, non-discrete, Hausdorff topological group.

  1. Any proper subgroup of Q/Z has an infinite index.
  2. A discrete subset of Q/Z cannot generate the entire group.

  From these two points, for any discrete subset A of Q/Z, the generated
  subgroup <A> is a proper subgroup, and thus the quotient group Q/Z / <A>
  has infinite cardinality.
  
  Therefore, the minimum possible cardinality I_{Q/Z} is countably infinite.
  Since the cardinality of the quotient group cannot exceed the cardinality of
  the group itself (which is countably infinite), this is the maximum possible value for I_G.
  
  The value is Aleph_0, or countably infinite. We represent this as "infinity".
  """
  
  # The largest value for I_G is countably infinite.
  largest_value = math.inf
  
  # We print a string representation as the answer is a cardinality.
  print("infinity")

solve()