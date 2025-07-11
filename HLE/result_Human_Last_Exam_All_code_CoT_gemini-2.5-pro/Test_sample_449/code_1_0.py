import math

# The starting position is (3000, 4000)
x = 3000
y = 4000

# Calculate the distance from the origin
r = math.sqrt(x**2 + y**2)

# The probability of the conditioned random walk starting at a large distance r
# from the origin ever hitting the set of immediate neighbors of the origin
# is approximated by 1 / log(r).
# The problem asks for the probability that it will *never* come to this set.
# This is 1 minus the hitting probability.

log_r = math.log(r)
prob_hit = 1 / log_r
prob_never_hit = 1 - prob_hit

# We need to output the equation with the numbers
print(f"The starting position is ({x}, {y}).")
print(f"The distance from the origin is r = sqrt({x}^2 + {y}^2) = {r}.")
print(f"The natural logarithm of the distance is log(r) = log({r}) = {log_r:.4f}.")
print(f"The probability of ever hitting the four neighbours of the origin is approximately 1 / log(r).")
print(f"P(hit) = 1 / {log_r:.4f} = {prob_hit:.4f}.")
print(f"The probability of never hitting the four neighbours is 1 - P(hit).")
print(f"P(never hit) = 1 - {prob_hit:.4f} = {prob_never_hit:.4f}.")
print(f"The result rounded to two significant digits is {prob_never_hit:.2f}.")
