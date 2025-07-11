import collections

# Step 1: Analyze the Data to find s_H and s_L, and e_H and e_L.
# The given probability distribution of (effort, output) pairs.
prob_dist = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375
}

# From y = s * e, we can find the possible values for the state 's'.
s_values = set()
for (e, y) in prob_dist.keys():
    s = y / e
    s_values.add(s)

s_list = sorted(list(s_values))
s_L = s_list[0]
s_H = s_list[1]

# Identify the two distinct effort levels from the data.
e_values = sorted(list(set(e for e, y in prob_dist.keys())), reverse=True)
e_H = e_values[0]
e_L = e_values[1]

print(f"Step 1: Inferred Model Parameters")
print(f"The two states of the world are s_L = {s_L} and s_H = {s_H}.")
print(f"The two effort levels chosen by the employee are e_L = {e_L} and e_H = {e_H}.\n")

# Step 2: Determine Employee's Beliefs (Conditional Expectations).
# The probability of choosing a certain effort is the marginal probability from the distribution.
# This corresponds to the probability of receiving the signal that induces that effort.
prob_e_H = prob_dist[(e_H, s_H * e_H)] + prob_dist[(e_H, s_L * e_H)]
prob_e_L = prob_dist[(e_L, s_H * e_L)] + prob_dist[(e_L, s_L * e_L)]

# P(s|e), which is equivalent to P(s|θ), can be found using the joint probabilities.
# P(s=s_H | e=e_H) = P(s=s_H, e=e_H) / P(e=e_H)
prob_sH_given_eH = prob_dist[(e_H, s_H * e_H)] / prob_e_H
prob_sL_given_eH = prob_dist[(e_H, s_L * e_H)] / prob_e_H

# P(s=s_H | e=e_L)
prob_sH_given_eL = prob_dist[(e_L, s_H * e_L)] / prob_e_L
prob_sL_given_eL = prob_dist[(e_L, s_L * e_L)] / prob_e_L

# The employee's conditional expectation of s for each effort/signal.
E_s_given_eH = s_H * prob_sH_given_eH + s_L * prob_sL_given_eH
E_s_given_eL = s_H * prob_sH_given_eL + s_L * prob_sL_given_eL

print(f"Step 2: Employee's Beliefs")
print(f"The employee's conditional expectation of s given the signal for high effort is E[s|θ_H] = {E_s_given_eH}.")
print(f"The employee's conditional expectation of s given the signal for low effort is E[s|θ_L] = {E_s_given_eL}.\n")


# Step 3: Calculate the Contract Parameter β.
# The employee's optimal effort is given by e = β * E[s|θ].
# We can solve for β using our data.
beta_from_eH = e_H / E_s_given_eH
beta_from_eL = e_L / E_s_given_eL

# In a consistent model, both calculations should yield the same β.
beta = beta_from_eH

print(f"Step 3: Calculating Contract Parameter β")
print(f"From the high effort choice, {e_H} = β * {E_s_given_eH}, we solve for β = {beta_from_eH:.1f}.")
print(f"From the low effort choice, {e_L} = β * {E_s_given_eL}, we solve for β = {beta_from_eL:.1f}.")
print(f"The value of the contract parameter β is consistently found to be {beta:.1f}.\n")

# Step 4: Solve the Firm's Problem for p.
# The firm chooses β to maximize its expected profit. The first-order condition
# from this optimization problem yields the relationship p = β.
p = beta

print(f"Step 4: Solving for Price p")
print("The firm's profit maximization problem implies that it will choose β such that p = β.")
print(f"The final equation is p = {beta:.1f}.")
print(f"Therefore, the value of p is {p:.1f}.")
