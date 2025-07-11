def solve():
  """
  The configuration space X_4 of 4-segment unit-length robot arms
  that form a closed loop is known to decompose into three disjoint
  connected manifolds without boundary.

  1. The main component has dimension 5. It can be viewed as a bundle
     of 3-spheres over a 2-sphere.
  2. The other two components are both diffeomorphic to the special
     orthogonal group SO(3), which is a 3-dimensional manifold.

  Therefore, the dimensions of the components are 5, 3, and 3.
  """
  y1 = 5
  y2 = 3
  y3 = 3
  
  # The question asks for the tuple of dimensions, sorted in descending order.
  # Although the problem uses ell for the number of components, which is 3,
  # the example format suggests printing the numbers separated by commas.
  print(f"{y1},{y2},{y3}")

solve()