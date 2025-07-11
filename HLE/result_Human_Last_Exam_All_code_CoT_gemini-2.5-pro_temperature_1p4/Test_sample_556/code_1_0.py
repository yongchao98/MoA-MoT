import collections

# Step 1: Define the provided data
# P(e, y) = probability
probabilities = {
    (22, 132): 0.4375,
    (22, 44): 0.0625,
    (10, 60): 0.0625,
    (10, 20): 0.4375,
}

# The two effort levels observed are e_H = 22 and e_L = 10.
# The employee chooses e_H=22 when receiving a high signal (theta_H) and e_L=10 for a low signal (theta_L).
e_H = 22
e_L = 10

# From y = s * e, we can determine the states s_H and s_L.
# Using the first data point:
s_H = 132 / e_H
# Using the second data point:
s_L = 44 / e_H

print(f"Derived high state of the world (s_H): {s_H}")
print(f"Derived low state of the world (s_L): {s_L}")
print("-" * 30)

# Step 2: Determine the contract parameter beta
# The employee's optimal effort is e = beta * E[s | theta].

# We can calculate E[s | theta_H] to find beta.
# First, find the probability of receiving the high signal, P(theta_H).
# This corresponds to all cases where e = 22.
P_theta_H = probabilities[(22, 132)] + probabilities[(22, 44)]

# Next, find the conditional probabilities of the states given the high signal.
# P(s_H | theta_H) = P(s_H and theta_H) / P(theta_H)
# P(s_H and theta_H) is the case where e=22 and s=s_H (y=132)
P_sH_given_thetaH = probabilities[(22, 132)] / P_theta_H
P_sL_given_thetaH = probabilities[(22, 44)] / P_theta_H

# Now, calculate the expected value of s given the high signal.
E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH

# Solve for beta using e_H = beta * E[s | theta_H]
beta = e_H / E_s_given_thetaH

print(f"The employee's optimal effort rule is e = beta * E[s|theta]")
print(f"Calculated conditional expectation E[s|theta_H] = {E_s_given_thetaH}")
print(f"Using e_H = {e_H}, the contract parameter beta is calculated to be: {beta}")
print("-" * 30)

# Step 3 & 4: Use the firm's optimality condition to find p.
# The firm chooses beta to maximize its expected profit.
# The expected profit Pi = E[p*y - w] = E[p*y - (alpha + beta*y)]
# Assuming the participation constraint is binding, the firm's problem simplifies to maximizing the total surplus:
# Pi(beta) = E[p*y(beta) - e(beta)^2/2]
# The first-order condition d(Pi)/d(beta) = 0 gives the optimal beta.
# After derivation (as shown in the plan), this FOC simplifies to a direct relationship:
# p * dE[y]/d(beta) = dE[e^2/2]/d(beta)
# which ultimately yields p = beta.

p = beta

print("The firm's profit maximization problem implies a direct relationship between price 'p' and the contract parameter 'beta'.")
print("The final equation is:")
print(f"p = beta")
print(f"Since beta = {beta}, the value of p must be {p}.")

# Final answer in the required format
final_answer = p
print(f"\nThe final calculated value of p is: {final_answer}")