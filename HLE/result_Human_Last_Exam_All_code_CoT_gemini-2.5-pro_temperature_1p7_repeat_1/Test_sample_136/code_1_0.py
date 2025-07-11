# In this problem, we need to find the limit of the conditional probability:
# P[x_0 was not visited | 0 was not visited]
# where the walk runs up to time t_n = n^2 * ln(n)^2.

# Let p be the probability that a random walk starting from x_0 returns to x_0 before hitting 0.
# P_{x_0}(T_{x_0}^+ < T_0)
# From the derivation, we found that x_0 has 4 neighbors.
# - 2 are common neighbors with 0 (c_1, c_2).
# - 2 are not (f_1, f_2).
# Let h(y) = P_y(T_{x_0} < T_0).
# The probability p is given by p = (1/4) * (h(c_1) + h(c_2) + h(f_1) + h(f_2)).

# By symmetry, for a walk starting at a common neighbor (on the bisector of 0 and x_0),
# the probability of hitting x_0 before 0 is 1/2.
h_c = 1/2

# For a walk starting at a neighbor f_i, which is very close to x_0 but far from 0,
# the probability of hitting x_0 first approaches 1 in the limit n -> infinity.
h_f = 1

# The number of common neighbors. From the problem description, this is 2.
k = 2

# The degree of any vertex in the 2D torus is 4.
d = 4

# Calculate p
p_numerator = k * h_c + (d - k) * h_f
p_denominator = d
p = p_numerator / p_denominator

# The final result is the square of this probability, p^2.
result_numerator = p_numerator**2
result_denominator = p_denominator**2
result_float = result_numerator / result_denominator

# Print the final equation
print(f"The limit is given by the square of a calculated escape probability.")
print(f"Let p be the probability that a walk from x_0 returns to x_0 before hitting 0.")
print(f"From the analysis of its neighbors, this probability is:")
print(f"p = (2 * (1/2) + 2 * 1) / 4 = ({int(2*h_c)} + {d-k}) / {d} = {int(p_numerator)} / {p_denominator}")
print(f"The limiting conditional probability is p^2.")
print(f"Limit = ({int(p_numerator)} / {p_denominator})^2 = {result_numerator} / {result_denominator} = {result_float}")