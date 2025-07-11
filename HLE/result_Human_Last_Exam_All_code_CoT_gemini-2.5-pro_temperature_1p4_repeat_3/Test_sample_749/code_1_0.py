import math

# Step 1: Define the hitting probabilities in the h->0 limit.
# As h->0, the environment becomes all-blue with probability approaching 1.
# In an all-blue environment, p_L = 1/5, p_R = 4/5.
# The hitting probability pi_k for a particle starting at k is:
# pi_k = 1 for k < 0 (drift is towards 0)
# pi_k = (p_L/p_R)^k = (1/4)^k for k > 0 (drift is away from 0)
pi_minus_1 = 1
pi_1 = 1/4

print(f"In the h->0 limit, we analyze the condition for infinite visits, which is R > 1.")
print(f"R is the reproduction number at site 0.")
print(f"R approaches (1+h) * (p_L(0)*pi_-1 + p_R(0)*pi_1).")
print(f"We have pi_-1 = {pi_minus_1} and pi_1 = {pi_1}.")
print("The jump probabilities p_L(0) and p_R(0) depend on the color of site 0.\n")

# Step 2: Analyze the case where site 0 is red.
# This occurs with probability h.
# If site 0 is red, p_L(0) = 4/5, p_R(0) = 1/5.
p_L_red = 4/5
p_R_red = 1/5
R_factor_red = p_L_red * pi_minus_1 + p_R_red * pi_1

print("--- Case 1: Site 0 is red (probability h) ---")
print(f"p_L(0) = {p_L_red}, p_R(0) = {p_R_red}")
print(f"The condition is (1+h) * ({p_L_red} * {pi_minus_1} + {p_R_red} * {pi_1}) > 1")
print(f"This simplifies to (1+h) * {R_factor_red} > 1")
h_threshold_red = 1/R_factor_red - 1
print(f"This requires 1+h > {1/R_factor_red:.4f}, so h > {h_threshold_red:.4f}.\n")

# Step 3: Analyze the case where site 0 is blue.
# This occurs with probability 1-h.
# If site 0 is blue, p_L(0) = 1/5, p_R(0) = 4/5.
p_L_blue = 1/5
p_R_blue = 4/5
R_factor_blue = p_L_blue * pi_minus_1 + p_R_blue * pi_1

print("--- Case 2: Site 0 is blue (probability 1-h) ---")
print(f"p_L(0) = {p_L_blue}, p_R(0) = {p_R_blue}")
print(f"The condition is (1+h) * ({p_L_blue} * {pi_minus_1} + {p_R_blue} * {pi_1}) > 1")
print(f"This simplifies to (1+h) * {R_factor_blue} > 1")
h_threshold_blue = 1/R_factor_blue - 1
print(f"This requires 1+h > {1/R_factor_blue:.4f}, so h > {h_threshold_blue:.4f}.\n")

# Step 4: Conclusion
print("--- Conclusion ---")
print(f"The probability of infinite visits is P(site 0 is red AND h>{h_threshold_red:.4f}) + P(site 0 is blue AND h>{h_threshold_blue:.4f}).")
print(f"This is h * I(h>{h_threshold_red:.4f}) + (1-h) * I(h>{h_threshold_blue:.4f}), where I is the indicator function.")
print("As we take the limit h -> 0, h is smaller than both thresholds.")
print("So the indicator functions are both 0.")
print("The probability converges to 0.")
