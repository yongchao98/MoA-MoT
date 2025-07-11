import collections

# Step 1: Deconstruct the Data and find the states of the world.
# The given data is P(e, y). The relationship is y = s * e, so s = y / e.
data = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375
}

# Calculate the states 's' for each (e, y) pair.
states = set()
for (e, y) in data.keys():
    s = y / e
    states.add(s)

# The states are s_H (high) and s_L (low).
s_H = max(states)
s_L = min(states)

print(f"Step 1: Identifying the states of the world from the data.")
print(f"The high state is s_H = {s_H}")
print(f"The low state is s_L = {s_L}\n")

# Step 2: Determine model probabilities.
# Let's associate higher effort (e=22) with the high signal (theta_H)
# and lower effort (e=10) with the low signal (theta_L).
e_H = 22
e_L = 10

# The probability distribution P(e,y) can be mapped to P(signal, state).
# P(theta_H, s_H) corresponds to P(e=22, y=132)
# P(theta_H, s_L) corresponds to P(e=22, y=44)
# P(theta_L, s_H) corresponds to P(e=10, y=60)
# P(theta_L, s_L) corresponds to P(e=10, y=20)
prob_joint = {
    ('theta_H', s_H): data[(e_H, s_H * e_H)],
    ('theta_H', s_L): data[(e_H, s_L * e_H)],
    ('theta_L', s_H): data[(e_L, s_H * e_L)],
    ('theta_L', s_L): data[(e_L, s_L * e_L)],
}

# Calculate marginal probabilities of signals
prob_theta_H = prob_joint[('theta_H', s_H)] + prob_joint[('theta_H', s_L)]
prob_theta_L = prob_joint[('theta_L', s_H)] + prob_joint[('theta_L', s_L)]

# Calculate conditional probabilities P(s | theta) using Bayes' rule: P(s|theta) = P(s, theta) / P(theta)
prob_sH_given_thetaH = prob_joint[('theta_H', s_H)] / prob_theta_H
prob_sL_given_thetaH = prob_joint[('theta_H', s_L)] / prob_theta_H
prob_sH_given_thetaL = prob_joint[('theta_L', s_H)] / prob_theta_L
prob_sL_given_thetaL = prob_joint[('theta_L', s_L)] / prob_theta_L

# Calculate the conditional expected value of the state, E[s | theta]
E_s_given_thetaH = s_H * prob_sH_given_thetaH + s_L * prob_sL_given_thetaH
E_s_given_thetaL = s_H * prob_sH_given_thetaL + s_L * prob_sL_given_thetaL

print("Step 2 & 3: Calculating conditional expectations and finding the bonus parameter 'β'.")
print(f"The expected state given the high signal is E[s | θ_H] = {s_H}*{prob_sH_given_thetaH:.3f} + {s_L}*{prob_sL_given_thetaH:.3f} = {E_s_given_thetaH}")
print(f"The expected state given the low signal is E[s | θ_L] = {s_H}*{prob_sH_given_thetaL:.3f} + {s_L}*{prob_sL_given_thetaL:.3f} = {E_s_given_thetaL}\n")

# Step 3: Find the bonus parameter beta from the employee's incentive compatibility (IC) condition: e = beta * E[s | theta]
# We can calculate beta from both high and low effort choices. They must be consistent.
beta_from_H = e_H / E_s_given_thetaH
beta_from_L = e_L / E_s_given_thetaL

# For the final result, we will just use one of them as they are identical.
beta = beta_from_H

print(f"The employee's optimal effort rule is e = β * E[s|θ].")
print(f"Using the high effort level: {e_H} = β * {E_s_given_thetaH}, which implies β = {beta_from_H:.2f}")
print(f"Using the low effort level: {e_L} = β * {E_s_given_thetaL}, which implies β = {beta_from_L:.2f}")
print(f"The determined value of the bonus parameter is β = {beta:.2f}\n")


# Step 4 & 5: Analyze the firm's problem to solve for p.
# The firm chooses beta to maximize profit: Π = p*E[y] - E[e^2/2].
# Let's express E[y] and E[e^2/2] in terms of beta.
# We know E[y] = E[s*e] = E[s*β*E[s|θ]] = β * E[s*E[s|θ]].
# and E[e^2/2] = E[(β*E[s|θ])^2 / 2] = (β^2 / 2) * E[(E[s|θ])^2].
# By the law of iterated expectations, E[s*E[s|θ]] = E[(E[s|θ])^2]. Let's call this constant K.
K = (E_s_given_thetaH**2) * prob_theta_H + (E_s_given_thetaL**2) * prob_theta_L
# So, Π(β) = p*β*K - (β^2/2)*K = K * (p*β - β^2/2)

# To maximize this profit, the firm takes the derivative with respect to beta and sets it to zero.
# dΠ/dβ = K * (p - β)
# Setting dΠ/dβ = 0 gives K*(p-β) = 0, which means p = β.

print("Step 4 & 5: Finding the price 'p' from the firm's profit maximization problem.")
print("The firm chooses β to maximize its profit. The first-order condition for profit maximization is:")
print("d(Profit)/d(β) = 0, which simplifies to the condition p = β.")
print(f"Since we determined from the data that β = {beta:.2f}, it must be that p = β.\n")

# Final result
p = beta
print("Final Calculation:")
print(f"p = β = {p:.2f}")
print("<<<4.0>>>")