import math

# The problem is about finding the probability that a random walk, starting at an angle of pi/2
# relative to the target axis, will eventually have its angle reach 0.
# The movement of the angle can be modeled as a 1D Brownian motion on the interval [0, 2*pi].
# This is a classic gambler's ruin problem.
# Let the starting position be x, on an interval [a, b].
# The probability of hitting 'a' before 'b' is (b-x)/(b-a).
# In our case, a = 0, b = 2*pi, and the starting angle x = pi/2.

a = 0
b = 2 * math.pi
x = math.pi / 2

# Probability of hitting 2*pi (the "upper" boundary) first is (x-a)/(b-a)
prob_hit_2pi = (x - a) / (b - a)

# Probability of hitting 0 (the "lower" boundary, which is the target direction) first
# is (b-x)/(b-a), or 1 - prob_hit_2pi.
prob_hit_0 = (b - x) / (b - a)

# Print the calculation steps
print(f"The starting angle is pi/2, which is {x} radians.")
print(f"The range of angles is [0, 2*pi], which is [{a}, {b}] radians.")
print("The probability of the angle reaching 0 before 2*pi is given by the formula (b-x)/(b-a).")
print(f"Calculation: ({b} - {x}) / ({b} - {a})")
print(f"= ({b-x}) / {b-a}")
print(f"= {1 - x/b}")

# The final probability calculation
final_prob_numerator = 3 * math.pi / 2
final_prob_denominator = 2 * math.pi
final_prob = final_prob_numerator / final_prob_denominator

# Print the equation representing the final answer
# The equation is 1 - (pi/2) / (2*pi) = 3/4
start_angle_str = "pi/2"
range_angle_str = "2*pi"
print(f"Final Equation: 1 - ({start_angle_str}) / ({range_angle_str}) = 1 - 1/4 = 3/4")
print(f"The final answer is {final_prob}")