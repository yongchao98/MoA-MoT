import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = captured_output = io.StringIO()

# --- Step 1: Deduce model parameters from the data ---
# Given probability distribution:
# P(e=22, y=132) = 0.4375
# P(e=22, y=44)  = 0.0625
# P(e=10, y=60)  = 0.0625
# P(e=10, y=20)  = 0.4375

# From y = s * e, we can find the states of the world, s.
# For e=22, y can be 132 or 44.
s1 = 132 / 22
s2 = 44 / 22
sH = max(s1, s2)
sL = min(s1, s2)
eH = 22
eL = 10

print(f"From the (e, y) pairs, we deduce the states of the world:")
print(f"High state s_H = {sH}")
print(f"Low state s_L = {sL}")
print("-" * 30)

# --- Step 2 & 3: Analyze employee behavior and calculate beta ---
# The employee chooses e = beta * E[s | theta].
# Effort e=22 corresponds to signal theta=s_H.
# Effort e=10 corresponds to signal theta=s_L.

# We map the given probabilities to P(theta, s)
P_thH_sH = 0.4375  # P(theta=s_H, s=s_H)
P_thH_sL = 0.0625  # P(theta=s_H, s=s_L)
P_thL_sH = 0.0625  # P(theta=s_L, s=s_H)
P_thL_sL = 0.4375  # P(theta=s_L, s=s_L)

# Calculate marginal probabilities for the signals
P_thH = P_thH_sH + P_thH_sL
P_thL = P_thL_sH + P_thL_sL

# Calculate conditional probabilities P(s | theta) using Bayes' rule
P_sH_given_thH = P_thH_sH / P_thH
P_sL_given_thH = P_thH_sL / P_thH
P_sH_given_thL = P_thL_sH / P_thL
P_sL_given_thL = P_thL_sL / P_thL

# Calculate conditional expectations E[s | theta]
E_s_given_thH = sH * P_sH_given_thH + sL * P_sL_given_thH
E_s_given_thL = sH * P_sH_given_thL + sL * P_sL_given_thL

print("The employee's optimal effort depends on their expectation of the state.")
print(f"Expected state given high signal E[s|theta=s_H] = {E_s_given_thH}")
print(f"Expected state given low signal E[s|theta=s_L] = {E_s_given_thL}")

# Solve for beta from e = beta * E[s | theta]
# We have two equations: eH = beta * E_s_given_thH and eL = beta * E_s_given_thL
beta_from_H = eH / E_s_given_thH
beta_from_L = eL / E_s_given_thL

# The value of beta must be consistent. We can use either one.
beta = beta_from_H

print(f"Using e_H = β * E[s|s_H], we get {eH} = β * {E_s_given_thH}, so β = {beta_from_H}")
print(f"Using e_L = β * E[s|s_L], we get {eL} = β * {E_s_given_thL}, so β = {beta_from_L}")
print(f"The contract parameter is β = {beta}")
print("-" * 30)

# --- Step 4 & 5: Analyze firm's behavior and solve for p ---
# The firm chooses beta to maximize Profit(beta) = E[p*y - w].
# A detailed derivation shows that the first-order condition for the firm's
# profit maximization problem simplifies to p = beta.

p = beta

print("The firm chooses β to maximize its profit.")
print("The first-order condition for this optimization yields the equation: p = β")
print(f"Since we found β = {beta}, the price p must be equal to this value.")
print(f"Final Equation: p = {beta}")
print("-" * 30)

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()
# Print the output to the original stdout
print(output)
print(f">>>{p}<<<")