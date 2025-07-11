# Given joint probabilities of (effort, output)
prob_e22_y132 = 0.4375
prob_e22_y44 = 0.0625
prob_e10_y60 = 0.0625
prob_e10_y20 = 0.4375

# From the data, we determine the states of the world (s = y/e) and effort levels.
# High state s_H = 132/22 = 60/10 = 6
# Low state s_L = 44/22 = 20/10 = 2
s_H = 6
s_L = 2

# High effort e_H corresponds to the signal for the high state.
e_H = 22

# The joint probabilities can be expressed in terms of (effort, state)
# P(e=22, s=6) = P(e=22, y=132)
prob_e22_sH = prob_e22_y132
# P(e=22, s=2) = P(e=22, y=44)
prob_e22_sL = prob_e22_y44

# To calculate E[s | e=22], we need the conditional probabilities P(s | e=22).
# First, calculate the marginal probability of observing the high effort, P(e=22).
prob_e22 = prob_e22_sH + prob_e22_sL

# Now, calculate the conditional probabilities using Bayes' rule: P(s|e) = P(s,e) / P(e)
prob_sH_given_e22 = prob_e22_sH / prob_e22
prob_sL_given_e22 = prob_e22_sL / prob_e22

# Calculate the expected value of the state, conditional on observing the high effort.
# This corresponds to E[s | θ = s_H].
E_s_given_e22 = s_H * prob_sH_given_e22 + s_L * prob_sL_given_e22

# From the employee's optimal effort rule, we have e_H = β * E[s | e=22].
# We can solve for β.
beta = e_H / E_s_given_e22

# As derived in the explanation, the firm's optimal choice leads to p = β.
p = beta

# Print the steps of the calculation for beta
print(f"The high effort level is e_H = {e_H}.")
print(f"The states of the world are s_H = {s_H} and s_L = {s_L}.")
print(f"The marginal probability of the high effort is P(e={e_H}) = {prob_e22_sH} + {prob_e22_sL} = {prob_e22}.")
print(f"The conditional probability of the high state given high effort is P(s={s_H}|e={e_H}) = {prob_e22_sH} / {prob_e22} = {prob_sH_given_e22:.4f}.")
print(f"The conditional probability of the low state given high effort is P(s={s_L}|e={e_H}) = {prob_e22_sL} / {prob_e22} = {prob_sL_given_e22:.4f}.")
print(f"The expected state given high effort is E[s|e={e_H}] = {s_H}*{prob_sH_given_e22:.4f} + {s_L}*{prob_sL_given_e22:.4f} = {E_s_given_e22}.")
print(f"From the employee's effort rule, β = e_H / E[s|e={e_H}].")
print(f"So, β = {e_H} / {E_s_given_e22} = {beta}.")
print(f"Since the firm chooses β optimally, we have p = β.")
print(f"Therefore, the value of p is {p}.")
