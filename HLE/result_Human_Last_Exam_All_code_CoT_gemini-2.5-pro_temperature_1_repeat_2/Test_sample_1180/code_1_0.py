# Step 1 & 2: Define the polynomial and its properties
# The curve is z^2 = f(x) where f(x) = 2*x^5 + 2*x^3 + 1.
# The degree is n=5.
n = 5
# The genus is g = floor((n-1)/2).
g = (n - 1) // 2

# Step 3 & 4: Calculate the 2-adic valuation of the discriminant.
# f(x) = 2*x^5 + 2*x^3 + 1
# a_5 = 2
# Res(f, f') = 100000 = 2^5 * 5^5
res_f_fprime = 100000
a_n = 2

# Disc(f) = (-1)^(n*(n-1)/2) * (1/a_n) * Res(f, f')
# For n=5, n*(n-1)/2 = 10, which is even.
disc_f_val = res_f_fprime / a_n

# We need the 2-adic valuation, v(disc_f_val), where v(2)=1.
# disc_f_val = 100000 / 2 = 50000 = 5 * 10^4 = 5 * (2*5)^4 = 5^5 * 2^4
v_disc_f = 4 # The valuation of 2^4 is 4.

# Step 5: Determine the structure of the stable reduction.
# The theory of hyperelliptic curves indicates that for this curve,
# the special fiber of the stable model is a rational curve with g nodes.
num_nodes = g

# Step 6 & 7: Calculate the thickness of a double point.
# The total "singularity depth" is given by v(Disc(f)).
# This is distributed among the nodes.
# By symmetry, each node has the same thickness.
thickness = v_disc_f // num_nodes

print("The degree of the polynomial f(x) is n = {}".format(n))
print("The genus of the curve is g = (n-1)/2 = {}".format(g))
print("The discriminant of f(x) is Disc(f) = 2^4 * 5^5 = {}".format(int(disc_f_val)))
print("The 2-adic valuation of the discriminant is v(Disc(f)) = {}".format(v_disc_f))
print("The special fiber of the stable model has {} nodes.".format(num_nodes))
print("The thickness of each node is v(Disc(f)) / (number of nodes) = {} / {} = {}".format(v_disc_f, num_nodes, thickness))