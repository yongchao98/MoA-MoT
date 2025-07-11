import math

# We need to find the minimum q (q_0) such that the probability of Alice winning (p_A) is at least 0.95.
# The relationship between p_A and q is given by the equation: p_A = 1 - (1 - q^4 * p_A^3)^3
# We set p_A = 0.95 and solve for q_0.

# The final equation to solve for q_0 is:
# q_0 = ((1 - 0.05**(1/3)) / (0.95**3))**(1/4)

# Let's calculate each part of the equation.
p_A = 0.95
win_prob_complement = 1 - p_A  # This is 0.05

# 1. Calculate the term for Bob forcing a loss: 0.05^(1/3)
# This represents the value of (1 - q_0^4 * p_A^3)
term1_cbrt = win_prob_complement**(1/3)
print(f"The value of (0.05)^(1/3) is: {term1_cbrt}")

# 2. Calculate the term related to Alice's probability at her turn: 0.95^3
# This is part of the expression for p_B = (q*p_A)^3
term2_pA_cubed = p_A**3
print(f"The value of (0.95)^3 is: {term2_pA_cubed}")

# 3. Calculate the numerator of the fraction: 1 - 0.05^(1/3)
# This corresponds to q_0^4 * p_A^3
numerator = 1 - term1_cbrt
print(f"The numerator (1 - 0.05^(1/3)) is: {numerator}")

# 4. Calculate q_0^4
# This is the ratio of the numerator and term2_pA_cubed
q0_pow_4 = numerator / term2_pA_cubed
print(f"The value of q_0^4 is: {q0_pow_4}")

# 5. Calculate q_0 by taking the fourth root
q0 = q0_pow_4**(1/4)
print(f"The minimum value of q (q_0) is: {q0}")

# 6. The problem asks for floor(100 * q_0)
final_answer = math.floor(100 * q0)

print("\nFinal calculation:")
print(f"100 * q_0 = {100 * q0}")
print(f"The floor of (100 * q_0) is: {final_answer}")

print(f"\n{final_answer}")