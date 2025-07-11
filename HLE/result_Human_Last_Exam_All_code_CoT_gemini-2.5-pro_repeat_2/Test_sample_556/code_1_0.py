import sys
import io

# Backup the original stdout
original_stdout = sys.stdout
# Create a new stream to capture the output
captured_output = io.StringIO()
# Redirect stdout to the new stream
sys.stdout = captured_output

# Step 1: Define variables from the problem description and data analysis.
# From the P(e,y) data, we can infer the states (s) and efforts (e).
# e=22, y=132 -> s = 132/22 = 6
# e=22, y=44  -> s = 44/22  = 2
# e=10, y=60  -> s = 60/10  = 6
# e=10, y=20  -> s = 20/10  = 2
s_H = 6.0
s_L = 2.0
# High effort is chosen with a high signal, low effort with a low signal.
e_H = 22.0 # Corresponds to signal θ=s_H
e_L = 10.0 # Corresponds to signal θ=s_L

# The problem gives the joint probabilities of outcomes. We map these to (s, θ) events.
# P(e=22, y=132) corresponds to P(s=s_H, θ=s_H)
P_sH_thetaH = 0.4375
# P(e=22, y=44) corresponds to P(s=s_L, θ=s_H)
P_sL_thetaH = 0.0625
# P(e=10, y=60) corresponds to P(s=s_H, θ=s_L)
P_sH_thetaL = 0.0625
# P(e=10, y=20) corresponds to P(s=s_L, θ=s_L)
P_sL_thetaL = 0.4375

# Step 2: Calculate marginal and conditional probabilities.
# Marginal probability of receiving the high signal θ=s_H
P_thetaH = P_sH_thetaH + P_sL_thetaH

# Conditional probabilities of the state, given the high signal
P_sH_given_thetaH = P_sH_thetaH / P_thetaH
P_sL_given_thetaH = P_sL_thetaH / P_thetaH

# Step 3: Calculate the employee's conditional expectation of the state.
# The employee's expectation of the state 's' when receiving the high signal θ=s_H
E_s_given_thetaH = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH

# Step 4: Calculate β from the employee's optimal effort rule: e_H = β * E[s|θ=s_H]
beta = e_H / E_s_given_thetaH

# Step 5: The firm's optimization implies p = β.
p = beta

# Step 6: Print the final calculation and the result.
# The final equation is p = β = e_H / E[s|θ=s_H]
print("The final price 'p' is derived from the firm's optimal choice of β.")
print("The firm's optimization results in p = β.")
print("β is determined by the employee's optimal effort choice: β = e_H / E[s | θ=s_H]")
print("\nCalculating the components:")
print(f"High effort e_H = {e_H}")
print(f"Expectation E[s | θ=s_H] = s_H * P(s=s_H|θ=s_H) + s_L * P(s=s_L|θ=s_H)")
print(f"                       = {s_H} * {P_sH_given_thetaH:.4f} + {s_L} * {P_sL_given_thetaH:.4f}")
print(f"                       = {s_H * P_sH_given_thetaH:.4f} + {s_L * P_sL_given_thetaH:.4f}")
print(f"                       = {E_s_given_thetaH}")
print(f"\nTherefore, β = {e_H} / {E_s_given_thetaH} = {beta}")
print(f"Since p = β, the value of p is: {p}")

# Restore the original stdout
sys.stdout = original_stdout
# Get the captured output as a string
output_string = captured_output.getvalue()
# Print the string to the real stdout
print(output_string)

# Final answer in the required format
print(f"<<<{p}>>>")