import sys
import io

# Backup original stdout
original_stdout = sys.stdout
# Create a new StringIO object
captured_output = io.StringIO()
# Redirect stdout
sys.stdout = captured_output


# --- Start of the solution code ---

# Define the joint probabilities from the problem statement
# P(e, y)
P_e22_y132 = 0.4375
P_e22_y44 = 0.0625
P_e10_y60 = 0.0625
P_e10_y20 = 0.4375

# Step 1: Determine states of the world (s_H, s_L) and effort levels (e_H, e_L)
print("Step 1: Determine the parameters from the data.")
# The two distinct effort levels are given in the data:
e_H = 22  # High effort
e_L = 10  # Low effort

# The output is given by y = s * e, so the state of the world is s = y / e.
# For e = 22, we find the possible values for s:
s_val1 = 132 / e_H
s_val2 = 44 / e_H
# For e = 10, we find the same possible values for s:
s_val3 = 60 / e_L
s_val4 = 20 / e_L

# The two states of the world are the unique values found, s_H > s_L
s_H = max(s_val1, s_val2)
s_L = min(s_val1, s_val2)
print(f"The effort levels are e_H = {e_H} and e_L = {e_L}.")
print(f"The states of the world are s_H = {s_H} and s_L = {s_L}.")
print("-" * 30)

# Step 2: Calculate conditional expectations E[s|θ]
print("Step 2: Calculate the conditional expectation of the state.")
# Let θ_H be the signal that leads to the high effort e_H = 22.
# Let θ_L be the signal that leads to the low effort e_L = 10.
# The given probabilities P(e, y) correspond to the joint probabilities P(θ, s).
P_thetaH_sH = P_e22_y132 # P(θ=s_H, s=s_H)
P_thetaH_sL = P_e22_y44  # P(θ=s_H, s=s_L)
P_thetaL_sH = P_e10_y60  # P(θ=s_L, s=s_H)
P_thetaL_sL = P_e10_y20  # P(θ=s_L, s=s_L)

# Calculate the marginal probability of the high signal, P(θ_H)
P_thetaH = P_thetaH_sH + P_thetaH_sL

# Calculate the conditional probabilities of the state, given the high signal
# P(s|θ_H) = P(s, θ_H) / P(θ_H)
P_sH_given_thetaH = P_thetaH_sH / P_thetaH
P_sL_given_thetaH = P_thetaH_sL / P_thetaH

# Calculate the conditional expectation E[s|θ_H]
E_s_given_theta_H = s_H * P_sH_given_thetaH + s_L * P_sL_given_thetaH
print(f"The marginal probability of the high signal P(θ_H) is {P_thetaH:.4f}.")
print(f"The conditional expectation of the state given the high signal, E[s|θ_H], is {E_s_given_theta_H:.2f}.")
print("-" * 30)


# Step 3: Calculate β
print("Step 3: Calculate the contract parameter β.")
# From the employee's optimization, we have the effort rule: e(θ) = β * E[s|θ].
# We can use the high-effort case to find β: e_H = β * E[s|θ_H]
beta = e_H / E_s_given_theta_H
print("The employee's effort rule is e = β * E[s|θ].")
print(f"Using the high effort case: β = e_H / E[s|θ_H]")
print(f"β = {e_H} / {E_s_given_theta_H:.2f} = {beta:.2f}")
print("-" * 30)


# Step 4: Determine p
print("Step 4: Determine the price p.")
# From the firm's profit maximization problem, the first-order condition simplifies to p = β.
p = beta
print("The firm's profit maximization leads to the condition p = β.")
print(f"Therefore, p = {p:.2f}.")
print("-" * 30)

# Final Answer
print("The final equation for p is:")
print(f"p = β = {e_H} / (s_H * P(s_H|θ_H) + s_L * P(s_L|θ_H))")
print(f"p = {e_H} / ({s_H} * {P_sH_given_thetaH:.4f} + {s_L} * {P_sL_given_thetaH:.4f})")
print(f"p = {e_H} / ({E_s_given_theta_H:.2f})")
print(f"p = {p:.2f}")

# --- End of the solution code ---

# Restore original stdout
sys.stdout = original_stdout
# Get the captured output
output = captured_output.getvalue()

# Print the output to the user
print(output)
final_answer = p
print(f"<<<{final_answer}>>>")