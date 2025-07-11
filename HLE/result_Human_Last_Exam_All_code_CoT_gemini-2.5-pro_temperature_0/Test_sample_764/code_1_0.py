import numpy as np
from sklearn.linear_model import lars_path
from sklearn.datasets import make_regression

# Plan:
# 1. Generate synthetic data where the number of features (p) is greater than
#    the number of samples (n), a condition where Lasso solutions may not be unique.
# 2. Compute the entire Lasso solution path. This gives us a sequence of
#    penalty values (λ, called 'alpha' in scikit-learn) and the corresponding
#    coefficient vectors (β).
# 3. For each coefficient vector on the path, we calculate its L1 norm (t = sum|β_j|).
# 4. Finally, we print an explanation of why the equivalence between the two Lasso
#    formulations is not strictly one-to-one, and display a table showing the
#    relationship between λ and t for several points along the computed path.

# 1. Generate synthetic data with p > n
n_samples, n_features = 50, 100
X, y = make_regression(n_samples=n_samples, n_features=n_features, n_informative=10, noise=20.0, random_state=42)

# 2. Compute the Lasso solution path
# lars_path returns alphas (λ values) and the corresponding coefficients.
# The returned alphas are in decreasing order.
alphas, _, coefs = lars_path(X, y, method='lasso', verbose=False)

# 3. Calculate the L1 norm (t) for each coefficient vector in the path
# The shape of coefs is (n_features, n_steps), so we sum over the first axis.
l1_norms = np.sum(np.abs(coefs), axis=0)

# 4. Print explanation and results
print("The statement that the two Lasso formulations are strictly equivalent is not always true.")
print("The equivalence relies on the mapping between the budget 't' and the penalty 'λ'.")
print("This mapping is not always one-to-one. Specifically, when the solution to the")
print("penalized problem is not unique (e.g., when features > samples), a single λ value")
print("can correspond to a range of solutions with different L1 norms (t).")
print("Because of this, the equivalence is considered 'true in general' but 'not always strictly true'.")
print("\nThe table below shows corresponding values of λ and t from a computed solution path:")
print("As λ increases, the L1 norm 't' of the coefficients is forced to decrease.")

print("\n{:<20} {:<20}".format("λ (Penalty)", "t (L1 Norm Budget)"))
print("-" * 40)

# Reverse the arrays to show λ increasing from (near) zero
reversed_alphas = alphas[::-1]
reversed_l1_norms = l1_norms[::-1]

# Print a few representative points from the path
# We select evenly spaced points from the computed path for demonstration
num_points_to_show = 6
indices = np.linspace(0, len(reversed_alphas) - 1, num_points_to_show, dtype=int)

for i in indices:
    print("{:<20.4f} {:<20.4f}".format(reversed_alphas[i], reversed_l1_norms[i]))
