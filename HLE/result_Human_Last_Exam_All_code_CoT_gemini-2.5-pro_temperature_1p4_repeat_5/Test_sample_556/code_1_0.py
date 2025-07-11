import collections

# Step 1: Define the given data and determine s_H and s_L.
# P(e, y) distribution
prob_dist = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375
}

# Efforts and Outputs
efforts = [22, 10]
outputs = [132, 44, 60, 20]

# Determine states s = y/e
states = set()
for (e, y) in prob_dist.keys():
    states.add(y / e)

s_H = max(states)
s_L = min(states)
e_H = max(efforts)
e_L = min(efforts)

print(f"Step 1: Determine states of the world")
print(f"The high state of the world is s_H = {s_H}")
print(f"The low state of the world is s_L = {s_L}")
print("-" * 30)

# Step 2: Map probabilities from P(e,y) to P(theta, s)
# P(theta=s_H, s=s_H) corresponds to P(e=e_H, y=s_H*e_H)
p_H_H = prob_dist[(e_H, s_H * e_H)]
# P(theta=s_H, s=s_L) corresponds to P(e=e_H, y=s_L*e_H)
p_H_L = prob_dist[(e_H, s_L * e_H)]
# P(theta=s_L, s=s_H) corresponds to P(e=e_L, y=s_H*e_L)
p_L_H = prob_dist[(e_L, s_H * e_L)]
# P(theta=s_L, s=s_L) corresponds to P(e=e_L, y=s_L*e_L)
p_L_L = prob_dist[(e_L, s_L * e_L)]

print("Step 2: Map probabilities")
print(f"P(signal=s_H, state=s_H) = {p_H_H}")
print(f"P(signal=s_H, state=s_L) = {p_H_L}")
print(f"P(signal=s_L, state=s_H) = {p_L_H}")
print(f"P(signal=s_L, state=s_L) = {p_L_L}")
print("-" * 30)

# Step 3: Determine the contract parameter beta
# Marginal probability of signals
p_theta_H = p_H_H + p_H_L
p_theta_L = p_L_H + p_L_L

# Conditional expectation of s given theta
E_s_given_theta_H = (s_H * p_H_H + s_L * p_H_L) / p_theta_H
E_s_given_theta_L = (s_H * p_L_H + s_L * p_L_L) / p_theta_L

# Calculate beta from e(theta) = beta * E[s|theta]
beta = e_H / E_s_given_theta_H

print("Step 3: Determine contract parameter beta")
print(f"The employee's optimal effort rule is e = beta * E[s|signal].")
print(f"E[s|signal=s_H] = {E_s_given_theta_H}")
print(f"E[s|signal=s_L] = {E_s_given_theta_L}")
print(f"From e_H = beta * E[s|signal=s_H] => {e_H} = beta * {E_s_given_theta_H}")
print(f"beta = {beta}")
# Verification
# print(f"Verification: {e_L} = {beta} * {E_s_given_theta_L} => {e_L} = {beta * E_s_given_theta_L}")
print("-" * 30)


# Step 4: Determine the contract parameter alpha
# Calculate E[y] and E[e^2/2]
E_y = sum(y * p for (e, y), p in prob_dist.items())
E_e_sq_half = ((e_H**2 / 2) * p_theta_H) + ((e_L**2 / 2) * p_theta_L)

# The participation constraint E[w - e^2/2] >= 0 is binding, so E[alpha + beta*y - e^2/2] = 0
# alpha + beta*E[y] - E[e^2/2] = 0
alpha = E_e_sq_half - beta * E_y

print("Step 4: Determine contract parameter alpha")
print("Assuming the participation constraint is binding (employee's expected utility = 0):")
print(f"alpha + beta * E[y] - E[e^2/2] = 0")
print(f"E[y] = {E_y}")
print(f"E[e^2/2] = {E_e_sq_half}")
print(f"alpha + {beta} * {E_y} - {E_e_sq_half} = 0")
print(f"alpha = {E_e_sq_half} - {beta * E_y} = {alpha}")
print("-" * 30)

# Step 5: Determine the price p
# The zero-profit condition is E[p*y - w] = 0
# E[(p-beta)*y - alpha] = 0
# (p-beta)*E[y] - alpha = 0
# p = beta + alpha / E[y]
p = beta + alpha / E_y

print("Step 5: Determine the price p")
print("Assuming firms compete and expected profit is zero:")
print("E[Profit] = (p - beta) * E[y] - alpha = 0")
print(f"(p - {beta}) * {E_y} - ({alpha}) = 0")
print(f"(p - {beta}) * {E_y} = {alpha}")
print(f"p - {beta} = {alpha / E_y}")
print(f"p = {beta} + ({alpha / E_y})")
print(f"p = {p}")
print("-" * 30)
print(f"The final value of p is: {p}")