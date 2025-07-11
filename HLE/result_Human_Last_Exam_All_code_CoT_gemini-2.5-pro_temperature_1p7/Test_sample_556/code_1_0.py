import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# Step 1: Input data from the problem
# P(e, y) from the problem description
P_e22_y132 = 0.4375
P_e22_y44 = 0.0625
P_e10_y60 = 0.0625
P_e10_y20 = 0.4375

# Effort levels observed in the data
e_H = 22.0
e_L = 10.0

print("Step 1: Parse the problem data.")
print(f"We are given four outcomes for (effort, output) and their probabilities:")
print(f"P(e={e_H}, y=132) = {P_e22_y132}")
print(f"P(e={e_H}, y=44) = {P_e22_y44}")
print(f"P(e={e_L}, y=60) = {P_e10_y60}")
print(f"P(e={e_L}, y=20) = {P_e10_y20}")
print("-" * 30)

# Step 2: Determine the states of the world, s_H and s_L
# The output y is determined by y = s * e, which implies s = y / e.
# We can find the two states by looking at the outcomes for a given effort, e.g., e_H.
s_H = 132 / e_H
s_L = 44 / e_H

print("Step 2: Determine the states of the world (s_H and s_L).")
print("The employee's output is given by y = s * e, where s is the state of the world.")
print("We observe a high effort e_H = 22 and a low effort e_L = 10.")
print("When effort is e_H = 22, we see two outputs, y=132 and y=44. This allows us to calculate the two possible states:")
print(f"High State (s_H) = y / e_H = 132 / {e_H} = {s_H}")
print(f"Low State (s_L) = y / e_H = 44 / {e_H} = {s_L}")
print("-" * 30)

# Step 3: Establish the joint probability distribution of signals (theta) and states (s)
# The high effort level e_H corresponds to the employee receiving a high signal (theta_H).
# The low effort level e_L corresponds to the employee receiving a low signal (theta_L).
P_thetaH_sH = P_e22_y132
P_thetaH_sL = P_e22_y44
P_thetaL_sH = P_e10_y60
P_thetaL_sL = P_e10_y20

print("Step 3: Map the data to a joint probability distribution P(signal, state).")
print("A higher effort (e=22) implies the employee received a 'high' signal (theta_H).")
print("A lower effort (e=10) implies the employee received a 'low' signal (theta_L).")
print("The state is revealed by the output for a given effort level.")
print(f"P(signal=theta_H, state=s_H) = {P_thetaH_sH}")
print(f"P(signal=theta_H, state=s_L) = {P_thetaH_sL}")
print(f"P(signal=theta_L, state=s_H) = {P_thetaL_sH}")
print(f"P(signal=theta_L, state=s_L) = {P_thetaL_sL}")
print("-" * 30)

# Step 4: Calculate the employee's conditional expectation of the state after receiving a high signal
# The total probability of receiving a high signal is the sum of joint probabilities over all states.
P_thetaH = P_thetaH_sH + P_thetaH_sL

# The conditional probability of a state given a signal is P(s | theta) = P(s, theta) / P(theta)
P_sH_given_thetaH = P_thetaH_sH / P_thetaH
P_sL_given_thetaH = P_thetaH_sL / P_thetaH

# The conditional expectation is the weighted average of the states, using conditional probabilities as weights.
E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH

print("Step 4: Calculate the employee's conditional expectation E[s | theta_H].")
print("First, find the total probability of receiving the high signal theta_H:")
print(f"P(theta_H) = P(theta_H, s_H) + P(theta_H, s_L) = {P_thetaH_sH} + {P_thetaH_sL} = {P_thetaH}")
print("Next, find the conditional probabilities of the state given the high signal:")
print(f"P(s=s_H | theta=theta_H) = P(theta_H, s_H) / P(theta_H) = {P_thetaH_sH} / {P_thetaH} = {P_sH_given_thetaH:.3f}")
print(f"P(s=s_L | theta=theta_H) = P(theta_H, s_L) / P(theta_H) = {P_thetaH_sL} / {P_thetaH} = {P_sL_given_thetaH:.3f}")
print("The employee's expected state value, given a high signal, is:")
print(f"E[s | theta_H] = s_H * P(s_H|theta_H) + s_L * P(s_L|theta_H) = {s_H} * {P_sH_given_thetaH:.3f} + {s_L} * {P_sL_given_thetaH:.3f} = {E_s_given_thetaH}")
print("-" * 30)

# Step 5: Determine the contract parameter beta using the employee's optimal effort rule
# The employee chooses effort 'e' to maximize E[w - e^2/2 | theta], which leads to the rule: e = beta * E[s|theta].
# Using the high effort case: e_H = beta * E[s | theta_H]
beta = e_H / E_s_given_thetaH

print("Step 5: Determine the contract parameter beta.")
print("The employee chooses effort 'e' to maximize their expected payoff, leading to the optimal effort rule: e = beta * E[s | signal].")
print("Using the case of the high signal and high effort, we solve for beta:")
print(f"e_H = beta * E[s | theta_H]")
print(f"{e_H} = beta * {E_s_given_thetaH}")
print(f"beta = {e_H} / {E_s_given_thetaH} = {beta}")
print("-" * 30)

# Step 6: Relate beta to the price p using the firm's profit maximization problem
# The firm chooses beta to maximize its profit. In this principal-agent framework, the profit-maximizing choice for the firm is to set beta = p.
# This aligns the employee's incentive (to maximize beta * y) with the firm's goal (to maximize p * y).
p = beta

print("Step 6: Determine the price p from beta.")
print("The firm chooses the contract parameters to maximize its own expected profit.")
print("For this type of economic model, the firm's profit-maximizing strategy is to set the commission rate 'beta' equal to the output price 'p'.")
print("This makes the employee the 'residual claimant' on the output's value, providing the strongest incentives.")
print("The final equation is:")
print(f"p = beta = {p}")
print("-" * 30)

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output
output_str = captured_output.getvalue()
print(output_str)

# Final answer in the required format
final_answer = p
print(f'<<<{final_answer}>>>')