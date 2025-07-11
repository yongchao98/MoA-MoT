import numpy as np

# --- Problem Parameters ---
d = 101
alpha = 3
beta = 2

print("Step 1: Define the target quantity")
print("Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2")
print("-" * 30)

print("Step 2: Calculate E[||v||^2]")
print("The vector v is a reflection of d, and d is constructed to be a unit vector (||d||=1).")
print("Therefore, ||v|| = ||d|| = 1 for all samples.")
E_norm_v_sq = 1
print(f"E[||v||^2] = {E_norm_v_sq}")
print("-" * 30)

print("Step 3: Calculate E[v]")
print("E[v] = H * E[d], where H is the Householder reflection matrix.")
print("We first need to compute E[d].")

# Calculate E[d_1]
# For a ~ gamma(alpha, theta) and b ~ gamma(beta, theta), the ratio x = a/(a+b) ~ Beta(alpha, beta)
E_x = alpha / (alpha + beta)
E_d1 = 2 * E_x - 1
print(f"The first component of d is d_1 = (a-b)/(a+b).")
print(f"E[a/(a+b)] for a~Gamma({alpha},θ) and b~Gamma({beta},θ) is the mean of Beta({alpha},{beta}), which is {alpha}/({alpha}+{beta}) = {E_x:.2f}.")
print(f"E[d_1] = 2 * E[a/(a+b)] - 1 = 2 * {E_x:.2f} - 1 = {E_d1:.2f}")

# Explain E[d_rest]
print("The other components of d involve the term c/||c||, where c ~ N(0, I).")
print("Due to the symmetry of the Normal distribution, E[c/||c||] = 0.")
print("Therefore, the expectation of all other components of d is 0.")
print(f"So, E[d] = [{E_d1:.2f}, 0, ..., 0]^T")
print()

print("Now, we compute E[v] = H * E[d].")
print("H depends on (v1-v2). v1=e1, v2=1_d, so v1-v2 = [0, -1, ..., -1]^T.")
print("E[d] is a vector proportional to e1, i.e., E[d] = E[d_1] * e1.")
print("The product H*e1 simplifies because (v1-v2) is orthogonal to e1, i.e., (v1-v2)^T * e1 = 0.")
print("This leads to H*e1 = e1.")
print(f"So, E[v] = H * (E[d_1]*e1) = E[d_1] * (H*e1) = E[d_1] * e1 = E[d].")
print(f"E[v] = [{E_d1:.2f}, 0, ..., 0]^T")
print("-" * 30)


print("Step 4: Final Calculation")
norm_E_v_sq = E_d1**2
print(f"||E[v]||^2 = {E_d1:.2f}^2 = {norm_E_v_sq:.4f}")

trace_cov_v = E_norm_v_sq - norm_E_v_sq
print("\nPutting it all together:")
print(f"Tr(Cov(v)) = E[||v||^2] - ||E[v]||^2")
# The final line below prints out the numbers in the final equation.
print(f"Tr(Cov(v)) = {E_norm_v_sq} - {norm_E_v_sq:.4f} = {trace_cov_v:.4f}")

print(f"\n<<<{trace_cov_v}>>>")