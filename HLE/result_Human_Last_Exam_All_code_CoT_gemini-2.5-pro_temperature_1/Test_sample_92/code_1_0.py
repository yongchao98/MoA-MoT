# The problem asks for the probability P that the marble escapes.
# Let P be the probability of escaping (reaching bin 2025) and (1-P) be the probability of melting (reaching bin 2024).

# The ratio of probabilities of escaping versus melting can be argued to be equal to the ratio of
# the probabilities of jumping directly to the escape bin versus the melt bin from any given position n < 2024.
# Let's consider the starting position n = 0.
# The distance to the escape bin is 2025 - 0 = 2025.
# The distance to the melt bin is 2024 - 0 = 2024.
# The probability of jumping by i is (1/3)^|i|.

# Probability of jumping directly to bin 2025 is (1/3)^2025.
# Probability of jumping directly to bin 2024 is (1/3)^2024.

# The ratio of these probabilities is ((1/3)^2025) / ((1/3)^2024) = 1/3.
# This ratio holds for any starting position n < 2024.

# Let's set up the equation for the final probability P.
# The ratio of the probability of escaping to the probability of melting is 1/3.
# P / (1 - P) = 1 / 3

# We can solve this equation for P.
# 3 * P = 1 * (1 - P)
# 3 * P = 1 - P
# 3 * P + P = 1
# 4 * P = 1
# P = 1 / 4

p_numerator = 1
p_denominator = 4

print("Let P be the probability of escaping.")
print("The ratio of the probability of escaping to the probability of melting is given by:")
print("P / (1 - P) = 1 / 3")
print("Solving for P:")
print("3 * P = 1 * (1 - P)")
print("3 * P = 1 - P")
print("4 * P = 1")
print(f"P = {p_numerator} / {p_denominator}")

final_prob = p_numerator / p_denominator
print(f"The final probability is: {final_prob}")
