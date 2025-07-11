import math

# Step 1: Define the problem parameters
# The starting point of the random walk
x0 = (3000, 4000)
# The set of neighbors of the origin is A.
# We need the probability that the walk never enters A.

# Step 2: Formulate the probability using the potential kernel a(x)
# The probability that the conditioned walk starting at x0 hits the set A of the
# origin's neighbours is given by P_hit = a(y)/a(x0), where y is any neighbor
# of the origin. The probability we want is P_never_hit = 1 - P_hit.

# Step 3: Define the values of the potential kernel
# For a simple random walk on Z^2, it is a known result that the potential
# kernel at a neighbor of the origin is 1.
a_neighbor = 1.0

# For a point x far from the origin, a(x) is approximated by the formula:
# a(x) ~ (2/pi) * log(||x||) + (2*gamma + log(8))/pi
# where gamma is the Euler-Mascheroni constant.

# Step 4: Calculate a(x0) for x0 = (3000, 4000)
norm_x0 = math.sqrt(x0[0]**2 + x0[1]**2)

# Constants
gamma = 0.5772156649
pi = math.pi

# The constant term in the asymptotic formula
constant_term = (2 * gamma + math.log(8)) / pi
# The value of a(x0)
a_x0 = (2 / pi) * math.log(norm_x0) + constant_term

# Step 5: Calculate the final probability and print the results
# P_never_hit = 1 - a_neighbor / a_x0
p_never_hit = 1 - a_neighbor / a_x0

print("The probability that the walk never comes to the set of the four neighbours of the origin is given by the equation:")
print("P(never hit) = 1 - a(neighbor) / a(start_point)")
print(f"               = 1 - {a_neighbor} / a({x0[0]},{x0[1]})")
print("")
print(f"The potential kernel at a neighbor of the origin is a(neighbor) = {a_neighbor:.4f}")
print(f"The potential kernel at the starting point is a({x0[0]},{x0[1]}) ≈ {a_x0:.4f}")
print("")
print("Final Calculation:")
print(f"P(never hit) ≈ 1 - {a_neighbor:.4f} / {a_x0:.4f}")
print(f"             ≈ {p_never_hit:.4f}")
print("")
print(f"The approximate answer with two significant digits is {p_never_hit:.2f}.")

# The required final answer format
final_answer_val = round(p_never_hit, 2)