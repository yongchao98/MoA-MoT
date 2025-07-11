import sys

# Step 1: Deconstruct the Provided Data and Find Model Parameters

# Given joint probabilities P(e, y)
P_e22_y132 = 0.4375
P_e22_y44 = 0.0625
P_e10_y60 = 0.0625
P_e10_y20 = 0.4375

# From the data, we can identify two effort levels. Let's call them e_H and e_L.
e_H = 22.0
e_L = 10.0

# The output is given by y = s * e. We can find the states of the world, s, from this.
# For e_H = 22:
s_from_y132 = 132 / e_H # 6.0
s_from_y44 = 44 / e_H   # 2.0
# For e_L = 10:
s_from_y60 = 60 / e_L   # 6.0
s_from_y20 = 20 / e_L   # 2.0

# The two states of the world are s_H and s_L.
s_H = 6.0
s_L = 2.0

print(f"Step 1: Analyzing the data")
print(f"The high state of the world is s_H = {s_H}")
print(f"The low state of the world is s_L = {s_L}\n")

# Calculate marginal probabilities of effort, which correspond to the probability of receiving each signal.
# Let's assume observing signal theta_H leads to effort e_H, and theta_L to e_L.
P_e_H = P_e22_y132 + P_e22_y44
P_e_L = P_e10_y60 + P_e10_y20

# Calculate the employee's conditional beliefs (posterior probabilities) about the state s given the signal.
# P(s | theta) = P(s | e) = P(e, y(s)) / P(e)
P_sH_given_thetaH = P_e22_y132 / P_e_H
P_sL_given_thetaH = P_e22_y44 / P_e_H
P_sH_given_thetaL = P_e10_y60 / P_e_L
P_sL_given_thetaL = P_e10_y20 / P_e_L

print(f"The employee's belief that the state is high after seeing a high signal is P(s=s_H|theta=s_H) = {P_sH_given_thetaH:.4f}")
print(f"The employee's belief that the state is low after seeing a high signal is P(s=s_L|theta=s_H) = {P_sL_given_thetaH:.4f}")
print(f"The employee's belief that the state is high after seeing a low signal is P(s=s_H|theta=s_L) = {P_sH_given_thetaL:.4f}")
print(f"The employee's belief that the state is low after seeing a low signal is P(s=s_L|theta=s_L) = {P_sL_given_thetaL:.4f}\n")


# Calculate the employee's expected value of the state, conditional on the signal.
E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
E_s_given_thetaL = s_H * P_sH_given_thetaL + s_L * P_sL_given_thetaL

print(f"The employee's expected value of s after a high signal is E[s|theta_H] = {E_s_given_thetaH}")
print(f"The employee's expected value of s after a low signal is E[s|theta_L] = {E_s_given_thetaL}\n")


# Step 2: Determine the contract parameter beta
# The employee's optimal effort rule is e(theta) = beta * E[s|theta].
# We can solve for beta using the observed effort and the calculated expected values.
beta_from_H = e_H / E_s_given_thetaH
beta_from_L = e_L / E_s_given_thetaL

# The value of beta should be consistent for both effort levels.
beta = beta_from_H

print(f"Step 2: Determining the contract parameter beta")
print(f"From the high effort case: beta = e_H / E[s|theta_H] = {e_H} / {E_s_given_thetaH} = {beta}")
print(f"From the low effort case: beta = e_L / E[s|theta_L] = {e_L} / {E_s_given_thetaL} = {beta_from_L}")
print(f"The value is consistent. The contract parameter is beta = {beta}\n")


# Step 3 & 4: Analyze the Firm's Problem and Solve for p
# The firm chooses beta to maximize its expected profit: E[Profit] = E[(p - beta) * y - alpha]
# The firm's optimization problem, after substituting the participation constraint, becomes maximizing:
# E[Profit] = A * p * beta - B * beta^2 - C
# where A, B, and C are constants derived from model parameters.
# The first-order condition for maximizing this with respect to beta is:
# d(E[Profit])/d(beta) = A * p - 2 * B * beta = 0
# A * p = 2 * B * beta
# In our specific model, it can be shown that this simplifies to p = beta.

# Let's show this derivation briefly:
# E[Profit] = (p - beta) * E[y] - alpha
# E[y] = E[s*e] = E[s * beta * E[s|theta]] = beta * E[s * E[s|theta]]
# The term E[s*E[s|theta]] is a constant with respect to beta. Let's call it K1. E[y] = K1 * beta
# The participation constraint is alpha + beta*E[y] - E[e^2/2] >= U_out
# E[e^2/2] = E[(beta*E[s|theta])^2 / 2] = (beta^2/2) * E[(E[s|theta])^2]. Let K2 = E[(E[s|theta])^2]/2.
# alpha = U_out - beta*E[y] + E[e^2/2] = U_out - K1*beta^2 + K2*beta^2 = U_out - (K1-K2)*beta^2
# Substitute back into profit:
# E[Profit] = (p-beta)*K1*beta - (U_out - (K1-K2)*beta^2)
# E[Profit] = p*K1*beta - K1*beta^2 - U_out + K1*beta^2 - K2*beta^2 = p*K1*beta - K2*beta^2 - U_out
# The derivative w.r.t beta is p*K1 - 2*K2*beta = 0 => p = (2*K2/K1)*beta
# The key insight is that for this specific utility function (c(e)=e^2/2), K1 = 2*K2, which results in p = beta.

print(f"Step 3 & 4: Solving for p")
print("The firm chooses beta to maximize its expected profit. Economic theory for this model shows that the firm's optimal choice of beta is equal to the price, p.")
print("The first-order condition of the firm's profit maximization problem simplifies to:")
print("p = beta")
print("\nFrom Step 2, we found that the value of beta is 4.")
p = beta
print(f"Therefore, the final equation is p = {int(beta)}")
print(f"The value of p is {int(p)}.")

sys.stdout.flush()
# The final answer is 4
print("<<<4>>>")