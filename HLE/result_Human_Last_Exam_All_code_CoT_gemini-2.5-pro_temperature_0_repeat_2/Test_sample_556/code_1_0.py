import sys
import io

# Redirect stdout to capture print output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

# --- Start of the solution code ---

# Step 0: Define the given data
# Joint probabilities P(e, y)
P_e22_y132 = 0.4375
P_e22_y44 = 0.0625
P_e10_y60 = 0.0625
P_e10_y20 = 0.4375

# Observed effort levels
e_H = 22
e_L = 10

# Step 1: Determine the states of the world, s_H and s_L
# The relationship is y = s * e. We can find s from the (e, y) pairs.
# For e = 22, y can be 132 or 44.
s1 = 132 / e_H
s2 = 44 / e_H
# For e = 10, y can be 60 or 20.
s3 = 60 / e_L
s4 = 20 / e_L

# The states are consistent, so we can define s_H and s_L.
s_H = max(s1, s2)
s_L = min(s1, s2)

print(f"From the data, we deduce the two states of the world:")
print(f"High state s_H = {s_H}")
print(f"Low state s_L = {s_L}")
print("-" * 40)

# Step 2: Calculate the conditional expectation of s for each effort level.
# The employee's effort choice e depends on the signal θ. We associate e_H with θ_H and e_L with θ_L.
# The employee's optimal effort is e = β * E[s|θ].
# We first need to calculate E[s|e_H] and E[s|e_L].

# Probability of each effort level
P_eH = P_e22_y132 + P_e22_y44
P_eL = P_e10_y60 + P_e10_y20

# Joint probabilities P(s, e)
P_sH_eH = P_e22_y132
P_sL_eH = P_e22_y44
P_sH_eL = P_e10_y60
P_sL_eL = P_e10_y20

# Conditional probabilities P(s|e) = P(s,e) / P(e)
P_sH_given_eH = P_sH_eH / P_eH
P_sL_given_eH = P_sL_eH / P_eH
P_sH_given_eL = P_sH_eL / P_eL
P_sL_given_eL = P_sL_eL / P_eL

# Conditional expectations E[s|e] = s_H * P(s_H|e) + s_L * P(s_L|e)
E_s_given_eH = s_H * P_sH_given_eH + s_L * P_sL_given_eH
E_s_given_eL = s_H * P_sH_given_eL + s_L * P_sL_given_eL

print(f"The employee's expected state, given they choose high effort (e={e_H}), is E[s|e_H] = {E_s_given_eH}")
print(f"The employee's expected state, given they choose low effort (e={e_L}), is E[s|e_L] = {E_s_given_eL}")
print("-" * 40)

# Step 3: Calculate the incentive parameter β from the effort choice rule e = β * E[s|e].
# Using the high effort case:
beta_from_eH = e_H / E_s_given_eH
print(f"From the high effort case, the equation is: {e_H} = β * {E_s_given_eH}")
print(f"Solving for β gives: β = {e_H} / {E_s_given_eH} = {beta_from_eH}")

# Using the low effort case to verify:
beta_from_eL = e_L / E_s_given_eL
print(f"\nFrom the low effort case, the equation is: {e_L} = β * {E_s_given_eL}")
print(f"Solving for β gives: β = {e_L} / {E_s_given_eL} = {beta_from_eL}")

# The value is consistent.
beta = beta_from_eH
print(f"\nThe value for β is consistently {beta}.")
print("-" * 40)

# Step 4: Relate β to the price p.
# The firm chooses β to maximize expected profit. For this model's structure, the
# optimal β for the firm is equal to the price p of the output.
# This is a standard result from principal-agent theory which ensures the agent
# is properly incentivized to maximize the total value created.
p = beta

print("In this principal-agent model, the firm's profit-maximizing strategy is to set the incentive parameter β equal to the output price p.")
print(f"Therefore, p = β.")
print(f"The value of p is {p}.")

# --- End of the solution code ---

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())

# Final answer
final_answer = p
print(f"<<<{final_answer}>>>")