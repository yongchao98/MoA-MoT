import numpy as np

def skew_symmetric(v):
    """
    Creates a skew-symmetric matrix from a 3-element vector.
    v: 3x1 numpy array.
    """
    return np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [v[1], v[0], 0]
    ])

# --- Problem Setup ---
# Let's assume a 6-DOF state (attitude, position) for this example.
# state_error = [e_delta, e_position]

# The attitude deviation vector 'delta' to be reset (e.g., from a filter update)
delta = np.array([0.08, -0.04, 0.02])

# The pre-reset covariance matrix Sigma (for a 6x6 state)
# We assume some initial uncertainty and correlations.
# Attitude part has higher uncertainty in yaw (z-axis)
# Position part has some uncertainty
# There is some cross-correlation between attitude and position errors
Sigma = np.diag([np.deg2rad(0.5)**2, np.deg2rad(0.5)**2, np.deg2rad(1.0)**2, # attitude variance
                 0.1**2, 0.1**2, 0.2**2])                                # position variance
# Add some cross-correlation for demonstration
Sigma[0, 4] = Sigma[4, 0] = -0.001
Sigma[1, 3] = Sigma[3, 1] = 0.001


# --- Covariance Reset Calculation ---

# The transformation is only applied to the attitude states and their cross-correlations.
# 1. Construct the Jacobian for the attitude error transformation
I_3x3 = np.identity(3)
hat_delta = skew_symmetric(delta)
# The Jacobian G_delta is based on the BCH approximation for the reset operation
G_delta = I_3x3 + 0.5 * hat_delta

# 2. Construct the full Jacobian G for the entire state vector
# It's block-diagonal because the position error definition is unchanged.
G = np.block([
    [G_delta,      np.zeros((3, 3))],
    [np.zeros((3, 3)), np.identity(3) ]
])


# 3. Compute the post-reset covariance matrix Sigma_prime
# The formula is Sigma' = G * Sigma * G^T
Sigma_prime = G @ Sigma @ G.T


# --- Output the results ---
print("--- Attitude Reset Covariance Transformation ---\n")
print(f"The attitude deviation being reset is δ = {delta}")
print("\nThis vector δ is used to form the Jacobian G_δ for the attitude block:")
print("G_δ = I + 1/2 * hat(δ)")
print("I = \n", I_3x3)
print("hat(δ) = \n", np.round(hat_delta, 4))
print("G_δ = \n", np.round(G_delta, 4))

print("\nThe covariance matrix Σ transforms according to Σ' = G * Σ * G^T, where G is the block-Jacobian.")
print(f"\nThe full transformation Jacobian G = \n", np.round(G, 4))


print("\n--- Final Equation and Result ---")
# Final equation formatting
# Extracting the attitude block of Sigma and Sigma_prime for clarity
Sigma_dd = Sigma[0:3, 0:3]
Sigma_do = Sigma[0:3, 3:6]
Sigma_oo = Sigma[3:6, 3:6]

Sigma_prime_dd = Sigma_prime[0:3, 0:3]
Sigma_prime_do = Sigma_prime[0:3, 3:6]
Sigma_prime_oo = Sigma_prime[3:6, 3:6]

print("The exact expression for the post-reset covariance Σ' is:")
print("\n[ Σ'_δδ  Σ'_δo ]   [ G_δ  0 ] [ Σ_δδ  Σ_δo ] [ G_δ^T  0   ]")
print(  "[          ] = [      ] [          ] [        ]")
print(  "[ Σ'_oδ  Σ'_oo ]   [  0   I ] [ Σ_oδ  Σ_oo ] [  0     I   ]")

print("\nWhere G_δ = (I + 1/2 * hat(δ))")

print("\n--- Numerical Example ---")
print("\nPre-reset Covariance Σ =\n", np.round(Sigma, 6))
print("\nPost-reset Covariance Σ' =\n", np.round(Sigma_prime, 6))

print("\nNotice the changes:")
print("- The attitude-attitude block Σ_δδ is updated.")
print("- The attitude-other cross-correlation block Σ_δo is updated.")
print("- The other-other block Σ_oo remains unchanged.")

<<<print(f"Σ'_δδ = (I + 0.5 * hat({delta})) * Σ_δδ * (I + 0.5 * hat({delta}))^T")>>>