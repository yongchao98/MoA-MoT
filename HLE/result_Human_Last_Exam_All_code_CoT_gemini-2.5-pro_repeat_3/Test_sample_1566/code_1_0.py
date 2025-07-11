# The problem is to determine the number of topologically distinct continua
# that satisfy two given properties. The solution relies on logical deduction
# from theorems in topology rather than numerical computation.

# Based on the analysis:
# 1. The definition of an "end point" implies the continuum is chainable.
# 2. A theorem states that chainable continua have at most two end points.
# 3. Property (1) states there is more than one end point.
# --> Conclusion: The number of end points must be exactly 2.
num_endpoints = 2

# 4. Property (2) states there are exactly two orbits.
# 5. This forces the two end points to be one orbit and the "interior"
#    to be the second orbit. This implies the interior is homogeneous.
# 6. A chainable continuum with two end points and a homogeneous interior
#    is topologically equivalent to a simple arc (e.g., [0, 1]).

# The number of such distinct topological types is therefore 1.
number_of_continua = 1

# The final equation is simply the statement of this result.
# We print each number involved in the final statement.
# Final statement: The number of continua is 1.
print(f"{number_of_continua}")
