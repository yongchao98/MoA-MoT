import math

# This script calculates the minimum value of q (q_0) for Alice to have a 95% win probability.
# The final formula we solve is q_0 = [(1 - (1 - x)^(1/3)) / x^3]^(1/4), where x = 0.95.

# The target probability for Alice to win is x.
x = 0.95
print(f"The equation for q_0 is derived from the target win probability x = {x}.")

# The equation uses several numerical values derived from x.
# We will print each of these values as required.
# The term (1-x) appears, which corresponds to the probability of Bob winning.
prob_bob_wins = 1 - x
print(f"The value of (1 - x) in the equation is: {prob_bob_wins:.2f}")

# The cube root of this term is taken.
cbrt_term = math.pow(prob_bob_wins, 1/3)
print(f"The value of (1 - x)^(1/3) is: {cbrt_term}")

# The numerator of the main fraction is 1 minus the cube root term.
numerator = 1 - cbrt_term
print(f"The numerator 1 - (1 - x)^(1/3) is: {numerator}")

# The denominator is x cubed.
denominator = math.pow(x, 3)
print(f"The denominator x^3 is: {denominator}")

# The term inside the 4th root is the ratio of the numerator and denominator. This is q_0^4.
q_pow_4 = numerator / denominator
print(f"The value of q_0^4 is: {q_pow_4}")

# q_0 is the 4th root of the above value.
q0 = math.pow(q_pow_4, 1/4)
print(f"The calculated value of q_0 is: {q0}")

# The problem asks for the floor of 100 times q_0.
value_to_floor = 100 * q0
print(f"The value of 100 * q_0 is: {value_to_floor}")

# The final result is the floor of this value.
final_answer = math.floor(value_to_floor)
print(f"\nThe final answer is floor(100 * q_0) = {final_answer}")