import collections

# Step 1: Define the given probability distribution and extract key parameters.
# P(e, y) = prob
prob_data = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375,
}

# Extract unique values of e and s from the data.
efforts = sorted(list(set(e for e, y in prob_data.keys())), reverse=True)
states = sorted(list(set(y / e for e, y in prob_data.keys())), reverse=True)

e_H, e_L = efforts
s_H, s_L = states

print(f"From the data, we identify:")
print(f"High effort e_H = {e_H}, Low effort e_L = {e_L}")
print(f"High state s_H = {s_H}, Low state s_L = {s_L}")
print("-" * 30)

# Step 2 & 3: Create the P(θ, s) table and calculate conditional expectations.
# The event e=e_H is equivalent to the event θ=s_H.
# The event e=e_L is equivalent to the event θ=s_L.
P_theta_s = {
    ('H', 'H'): prob_data[(e_H, s_H * e_H)],  # P(θ_H, s_H)
    ('H', 'L'): prob_data[(e_H, s_L * e_H)],  # P(θ_H, s_L)
    ('L', 'H'): prob_data[(e_L, s_H * e_L)],  # P(θ_L, s_H)
    ('L', 'L'): prob_data[(e_L, s_L * e_L)],  # P(θ_L, s_L)
}

# Calculate marginal probabilities of signals P(θ).
P_theta_H = P_theta_s[('H', 'H')] + P_theta_s[('H', 'L')]
P_theta_L = P_theta_s[('L', 'H')] + P_theta_s[('L', 'L')]

# Calculate conditional expectations E[s | θ].
E_s_given_theta_H = (s_H * P_theta_s[('H', 'H')] + s_L * P_theta_s[('H', 'L')]) / P_theta_H
E_s_given_theta_L = (s_H * P_theta_s[('L', 'H')] + s_L * P_theta_s[('L', 'L')]) / P_theta_L

print(f"The employee's optimal effort rule is e(θ) = β * E[s|θ].")
print(f"The conditional expectation of s given the high signal is E[s|θ_H] = {E_s_given_theta_H}")
print(f"The conditional expectation of s given the low signal is E[s|θ_L] = {E_s_given_theta_L}")
print("-" * 30)

# Step 4: Calculate β from the employee's behavior.
# We use the equation for the high effort case: e_H = β * E[s | θ_H]
beta = e_H / E_s_given_theta_H

print("We solve for β using the high effort case:")
print(f"Equation: e_H = β * E[s | θ_H]")
print(f"Substituting values: {e_H} = β * {E_s_given_theta_H}")
print(f"β = {e_H} / {E_s_given_theta_H} = {beta}")

# As a check, we can use the low effort case.
# beta_check = e_L / E_s_given_theta_L
# print(f"Check with low effort case: β = {e_L} / {E_s_given_theta_L} = {beta_check}")
print("-" * 30)

# Step 5 & 6: Use the firm's optimality condition to find p.
# For a risk-neutral agent and linear contract, the firm's optimal choice is β = p.
p = beta

print("The firm chooses β to maximize profit.")
print("With a risk-neutral agent and linear contract, the optimal choice for the firm is to set β equal to p.")
print("Therefore, p = β.")
print(f"The final value of p is {p}.")
