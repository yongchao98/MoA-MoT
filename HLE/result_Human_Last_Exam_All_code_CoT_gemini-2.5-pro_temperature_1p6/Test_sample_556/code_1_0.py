import sys
import io

# Capture original stdout
original_stdout = sys.stdout
# Create a string buffer
output_buffer = io.StringIO()
# Redirect stdout to the buffer
sys.stdout = output_buffer

# --- Main logic ---

# The plan is as follows:
# 1. Deduce the high and low states of the world (s_H, s_L) from the provided (e, y) data.
# 2. Use the provided probability distribution P(e,y) to calculate the conditional expectation
#    of the state given the signal that leads to high effort (e=22).
# 3. Use the employee's optimal effort rule, e = β * E[s|θ], to solve for the contract parameter β.
# 4. Use the firm's profit maximization problem to show that the optimal contract sets p = β.
# 5. Conclude the value of p and print the derivation.

# Step 0: Define the given data from the problem
# P(e,y)
prob_e22_y132 = 0.4375
prob_e22_y44 = 0.0625
# Effort and output pairs
e_high = 22
y_high_sh = 132  # Output when e=22 and state is high
y_high_sl = 44   # Output when e=22 and state is low

# Step 1: Deduce the states of the world, s_H and s_L
# The relationship is y = s * e, so s = y / e.
# High state (s_H) corresponds to the higher output for a given effort.
s_H = y_high_sh / e_high
# Low state (s_L) corresponds to the lower output for a given effort.
s_L = y_high_sl / e_high

# Step 2: Calculate the conditional expectation E[s|θ_H]
# Let θ_H be the signal that induces the high effort e_high = 22.
# The joint probabilities involving this signal are given:
# P(θ=θ_H, s=s_H) is equivalent to P(e=22, y=132)
# P(θ=θ_H, s=s_L) is equivalent to P(e=22, y=44)
# The marginal probability of receiving signal θ_H is the sum:
prob_theta_H = prob_e22_y132 + prob_e22_y44
# Using the definition of conditional probability, P(A|B) = P(A,B)/P(B):
prob_sH_given_thetaH = prob_e22_y132 / prob_theta_H
prob_sL_given_thetaH = prob_e22_y44 / prob_theta_H
# The conditional expectation of s given the signal θ_H is:
E_s_given_thetaH = s_H * prob_sH_given_thetaH + s_L * prob_sL_given_thetaH

# Step 3: Solve for the contract parameter β
# The employee's optimal effort is given by the rule e = β * E[s|θ].
# For the high effort level, we have e_high = β * E_s_given_thetaH.
# We can solve this equation for β.
beta = e_high / E_s_given_thetaH

# Step 4: Relate the price p to the parameter β
# The firm chooses β to maximize its expected profit. It can be shown that the
# first-order condition of the firm's optimization problem simplifies to p = β.

# Step 5: Conclude the value of p and print the derivation
p = beta

print("The price 'p' is determined by solving the principal-agent model.")
print("The key steps are:")
print("\n1. Determine the states of the world, s_H and s_L, from the data.")
print(f"   From y = s * e, we find:")
print(f"   s_H = {y_high_sh} / {e_high} = {s_H}")
print(f"   s_L = {y_high_sl} / {e_high} = {s_L}")

print("\n2. Calculate the employee's conditional expectation of the state, E[s|θ_H], where θ_H is the signal leading to e=22.")
print(f"   E[s|θ_H] = s_H * P(s_H|θ_H) + s_L * P(s_L|θ_H)")
print(f"            = {s_H:.1f} * ({prob_e22_y132} / {prob_theta_H}) + {s_L:.1f} * ({prob_e22_y44} / {prob_theta_H})")
print(f"            = {s_H:.1f} * {prob_sH_given_thetaH:.4f} + {s_L:.1f} * {prob_sL_given_thetaH:.4f}")
print(f"            = {s_H * prob_sH_given_thetaH:.4f} + {s_L * prob_sL_given_thetaH:.4f} = {E_s_given_thetaH:.1f}")


print("\n3. Solve for the contract parameter β using the employee's effort rule: e = β * E[s|θ].")
print(f"   {e_high} = β * {E_s_given_thetaH:.1f}")
print(f"   β = {e_high} / {E_s_given_thetaH:.1f} = {beta:.1f}")


print("\n4. The firm's profit maximization implies an optimal contract where p = β.")
print("   Therefore, the value of p is:")
print(f"   p = {p:.1f}")

# --- End of Main logic ---

# Restore original stdout
sys.stdout = original_stdout
# Get the buffer's content
output_str = output_buffer.getvalue()

# Print the content to the real stdout
print(output_str)

final_answer = p
print(f'<<<{final_answer}>>>')