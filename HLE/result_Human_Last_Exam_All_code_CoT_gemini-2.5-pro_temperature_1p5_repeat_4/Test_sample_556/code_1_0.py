import collections

# Step 1: Define the given data and deduce the states of the world.
print("Step 1: Deducing states of the world from the data.")

# The given joint probability distribution of (effort, output) pairs.
P_ey = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375,
}

# The production function is y = s * e. We can find the states 's' from this.
s_values = set()
for (e, y) in P_ey.keys():
    # Calculate s for each (e, y) pair.
    s = y / e
    s_values.add(s)

# Identify the high and low states.
s_H = max(s_values)
s_L = min(s_values)
print(f"The states of the world are: s_H = {s_H} and s_L = {s_L}\n")

# Step 2: Analyze agent behavior and calculate probabilities.
print("Step 2: Calculating key probabilities from the data.")

# In the model, higher effort corresponds to the signal predicting a higher state.
e_H = 22  # Effort when signal is high (theta_H)
e_L = 10  # Effort when signal is low (theta_L)

# We can now map the joint distribution P(e, y) to P(e, s).
P_es = {
    (e_H, s_H): P_ey[(e_H, s_H * e_H)],
    (e_H, s_L): P_ey[(e_H, s_L * e_H)],
    (e_L, s_H): P_ey[(e_L, s_H * e_L)],
    (e_L, s_L): P_ey[(e_L, s_L * e_L)],
}

# The agent's choice of effort depends on the signal (theta).
# So, P(e=e_H) is the probability of receiving the high signal, P(theta=s_H).
P_theta_H = P_es[(e_H, s_H)] + P_es[(e_H, s_L)]
P_theta_L = P_es[(e_L, s_H)] + P_es[(e_L, s_L)]
print(f"Probability of high signal P(theta=s_H) = {P_theta_H}")
print(f"Probability of low signal P(theta=s_L) = {P_theta_L}\n")


# Step 3: Calculate the agent's conditional expectations.
print("Step 3: Calculating the agent's conditional expectations.")
# We need P(s|theta), which can be found using P(s|theta) = P(s, theta) / P(theta).
# Note that P(s, theta) is equivalent to P(s, e) in this context.
P_sH_given_thetaH = P_es[(e_H, s_H)] / P_theta_H
P_sL_given_thetaH = P_es[(e_H, s_L)] / P_theta_H

P_sH_given_thetaL = P_es[(e_L, s_H)] / P_theta_L
P_sL_given_thetaL = P_es[(e_L, s_L)] / P_theta_L

# Now, calculate the conditional expectation of the state, given the signal.
E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

print(f"Agent's expectation of s given a high signal, E[s|theta_H] = {E_s_given_thetaH}")
print(f"Agent's expectation of s given a low signal, E[s|theta_L] = {E_s_given_thetaL}\n")


# Step 4: Solve for the contract parameter beta.
print("Step 4: Solving for the incentive parameter beta.")
# The agent's optimal effort is e = beta * E[s|theta].
# We can find beta from both high and low effort choices.
beta_from_eH = e_H / E_s_given_thetaH
beta_from_eL = e_L / E_s_given_thetaL

# Both calculations should yield the same beta.
beta = beta_from_eH
print(f"From e_H = beta * E[s|theta_H] => {e_H} = beta * {E_s_given_thetaH} => beta = {beta_from_eH}")
print(f"From e_L = beta * E[s|theta_L] => {e_L} = beta * {E_s_given_thetaL} => beta = {beta_from_eL}")
print(f"The incentive parameter is beta = {beta}\n")


# Step 5 & 6: Analyze the firm's problem to find the price p.
print("Step 5 & 6: Solving for the price p from the firm's optimization.")
# The firm chooses beta to maximize its expected profit.
# The profit maximization problem simplifies, leading to the first-order condition:
# d(Profit)/d(beta) = p * K - beta * K = 0, where K is a term that cancels out.
# This gives the simple relationship: p = beta.

# The reasoning is that the firm's profit function is of the form:
# Profit(beta) = p*A*beta - 0.5*A*beta^2 - C
# where A = E[(E[s|theta])^2] and C is the agent's outside option utility.
# The derivative with respect to beta is p*A - A*beta. Setting this to 0 gives p = beta.

print("The firm's profit maximization implies a direct relationship between p and beta.")
p = beta
print(f"The final equation is p = beta.")
print(f"Since we calculated beta = {beta}, the price p must be {p}.")

print("\nFinal Answer:")
# The final result is the value of p
print(p)
>>> 4
