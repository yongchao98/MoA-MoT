import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import lasso_path
from sklearn.datasets import make_regression
from sklearn.preprocessing import StandardScaler

# Generate synthetic data for regression
# Use a well-posed problem (N > p) to demonstrate the general case
X, y, true_coef = make_regression(n_samples=100, n_features=10, n_informative=5, noise=20, coef=True, random_state=42)

# It's standard practice for Lasso to work with standardized variables
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# The `lasso_path` function computes the Lasso solution for a sequence of alpha values.
# In scikit-learn's terminology, `alpha` corresponds to `λ` in the problem statement.
alphas, coefs, _ = lasso_path(X_scaled, y, n_alphas=200, eps=1e-4)

# The `coefs` array has shape (n_features, n_alphas).
# We calculate the L1 norm for each set of coefficients along the path.
# This corresponds to 't' in the constrained problem formulation.
t_values = np.sum(np.abs(coefs), axis=0)

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# --- Plot 1: Coefficient paths vs. log(λ) ---
# We use a log scale for lambda (alpha) for better visualization of the path.
log_alphas = np.log10(alphas)
colors = plt.cm.viridis(np.linspace(0, 1, X.shape[1]))

for i in range(coefs.shape[0]):
    ax1.plot(log_alphas, coefs[i, :], color=colors[i], label=f'β_{i+1}')

ax1.set_xlabel('log(λ) [Penalty Strength]')
ax1.set_ylabel('Coefficient Value')
ax1.set_title('Lasso Coefficient Paths vs. Penalty (λ)')
ax1.grid(True)
ax1.legend(loc='upper right', title='Coefficients')
ax1.annotate('Path is continuous', xy=(-2, 5), xytext=(-3, 20),
             arrowprops=dict(facecolor='black', shrink=0.05),
             ha='center')


# --- Plot 2: L1 Norm (t) vs. log(λ) ---
ax2.plot(log_alphas, t_values, marker='.', linestyle='-', color='crimson')
ax2.set_xlabel('log(λ) [Penalty Strength]')
ax2.set_ylabel('L1 Norm of Coefficients (t = Σ|β|)')
ax2.set_title('Constraint (t) vs. Penalty (λ)')
ax2.grid(True)
# By convention, paths are often shown with the penalty decreasing (L1 norm increasing)
ax2.invert_xaxis()
ax2.annotate('Continuous & Monotonic', xy=(-3, 100), xytext=(-4, 150),
             arrowprops=dict(facecolor='black', shrink=0.05),
             ha='center')


plt.suptitle('Demonstration of Lasso Solution Path Properties', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
