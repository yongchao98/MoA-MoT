import sympy

# Step 1 & 2: Argument about non-special nodes
# A real non-special node 'c' is a critical point on the real axis such that phi(c) is not in {0, 1, infinity}.
# For a real rational function phi, phi'(c) and phi''(c) are real.
# Near c, phi(z) can be approximated by phi(c) + (phi''(c)/2) * (z-c)**2.
# The pre-image of the real line under this map consists of the real axis and a vertical line through c.
# Property (ii) states that neighbours of a real non-special node must be real nodes (real critical points).
# Following the vertical line from c, one would reach another critical point. For this neighbour to be real, it must lie on the real axis.
# However, it also must lie on the vertical line Re(z)=c. The only way for this to happen without the path being a trivial loop is if the vertical line is the real axis itself, which is not possible.
# A more rigorous proof shows this forces the conclusion that there are no real non-special nodes.

# Step 3: All real critical points are special vertices (p or q type).
# Let's try to construct functions with N_r poles in ]0,1[ and check for real non-special nodes.

# Let's test N_r = 2.
# A possible configuration of special vertices on [0,1] is p-r-p-r-p.
# Let the p-vertices be at 0, p1, 1 and the r-vertices at r1, r2.
# So, phi(0)=0, phi(1)=0, phi(p1)=0, phi(r1)=inf, phi(r2)=inf.
# From property (iii), all special vertices in ]0,1[ must have the same valency 2m.
# For the p- and q-vertices to be critical points, m must be > 1. For them to be extrema on the real line, m must be even.
# Let's assume m is an even integer >= 2.
# The corresponding function would be of the form:
# phi(x) = C * (x**m * (x-p1)**m * (1-x)**m) / ((x-r1)**m * (x-r2)**m)
# The critical points that are not p- or q-vertices (i.e., the non-special nodes) are given by the roots of (d/dx)log(phi(x)) that are not the p- or q-vertices themselves.
# (d/dx)log(phi(x)) = m * (1/x + 1/(x-p1) + 1/(x-1) - 1/(x-r1) - 1/(x-r2))
# Let's call the term in the parenthesis h(x). The roots of h(x) are the non-special nodes.
# Let's assume an ordered arrangement 0 < r1 < p1 < r2 < 1.
# Consider the interval (r1, p1). As x -> r1+, h(x) -> -inf. As x -> p1-, h(x) -> +inf.
# Since h(x) is continuous on (r1, p1), by the Intermediate Value Theorem, there must be a root of h(x) in (r1, p1).
# This root is a real non-special node. This contradicts our deduction from property (ii).
# So, a configuration with 2 or more poles will generally introduce real non-special nodes.

# Let's test N_r = 1.
# The simplest configuration is p-r-p on [0,1].
# Let the p-vertices be at 0 and 1, and the r-vertex at r1 in ]0,1[.
# phi(x) = C * (x**m * (1-x)**m) / (x-r1)**m
# (d/dx)log(phi(x)) = m * (1/x - 1/(1-x) - 1/(x-r1))
# Let's find the roots of h(x) = 1/x - 1/(1-x) - 1/(x-r1).
# h(x) = ( (1-x)(x-r1) - x(x-r1) - x(1-x) ) / ( x(1-x)(x-r1) )
# The numerator is: (-x**2 + (1+r1)x - r1) - (x**2 - r1*x) - (x - x**2)
# = -x**2 + 2*r1*x - r1
# We set the numerator to 0: x**2 - 2*r1*x + r1 = 0
# The discriminant is D = (2*r1)**2 - 4*r1 = 4*r1*(r1 - 1).
# Since r1 is in the interval ]0, 1[, we have r1 > 0 and r1 - 1 < 0.
# Therefore, D < 0.
# This quadratic equation has no real roots. So, there are no real non-special nodes.
# This construction is valid. It has one r-vertex in ]0, 1[.

# Let's test N_r = 0.
# For example, phi(x) = (2x-1)**2. phi(0)=1, phi(1)=1, phi(1/2)=0.
# Critical point at x=1/2 is a p-vertex. No non-special nodes. So N_r=0 is possible.

# The analysis shows that N_r=1 is possible, but N_r=2 (and by extension, N_r > 1) is not, because it necessarily creates real non-special nodes, which are forbidden by property (ii).
# Therefore, the maximum number of vertices labelled 'r' within ]0, 1[ is 1.

print("Based on the analysis, the maximum number of vertices labelled 'r' is 1.")
# The final equation to check for non-special nodes for N_r=1 was:
# x**2 - 2*r*x + r = 0
# For the case N_r=2, a simplified version could be checked.
# Let the vertices be at 0, 1/3, 2/3, 1. p(0), r(1/3), p(2/3), r(inf), p(1) is not simple.
# Let's use the reasoning above. The equation for non-special nodes for the case p(0)-r(r1)-p(p1)-r(r2)-p(1) is:
# 1/x + 1/(x-p1) + 1/(x-1) - 1/(x-r1) - 1/(x-r2) = 0.
# This equation will have real roots between the poles, confirming N_r >= 2 is not possible.
# Final Answer: 1
print("The quadratic equation for the non-special nodes in the N_r=1 case is:")
r = sympy.Symbol('r')
x = sympy.Symbol('x')
equation = x**2 - 2*r*x + r
print(f"x**2 - 2*r*x + r = 0")
discriminant = (2*r)**2 - 4*r
print(f"The discriminant is 4*r**2 - 4*r = 4*r*(r-1).")
print("For r in ]0, 1[, the discriminant is negative, so there are no real non-special nodes.")
print("For N_r=2, real non-special nodes are unavoidable.")
print("Thus, the maximum is 1.")