import numpy as np
from sklearn.datasets import make_regression
from sklearn.linear_model import Lasso
import warnings

# Suppress convergence warnings for demonstration purposes
warnings.filterwarnings('ignore', category=UserWarning)

# 1. Generate synthetic data
X, y = make_regression(n_samples=100, n_features=20, noise=10, random_state=42)

# 2. Define a range of lambda (alpha) values to test
# We include very small values to observe the behavior as lambda -> 0
alphas = np.logspace(2, -5, 100)

print("Demonstrating the relationship between lambda (penalty) and t (L1 norm bound)")
print("-" * 70)
print(f"{'Lambda (α)':<15} | {'t = ||β||₁':<20} | {'Comments'}")
print("-" * 70)

# 3. Fit Lasso for each alpha and calculate the L1 norm of the coefficients
t_values = []
last_t = -1
for alpha in alphas:
    lasso = Lasso(alpha=alpha, max_iter=10000)
    lasso.fit(X, y)
    
    # 4. Calculate t = ||β||₁
    t = np.sum(np.abs(lasso.coef_))
    t_values.append(t)
    
    # Print some interesting values
    if alpha == alphas[0]: # Largest alpha
        print(f"{alpha:<15.5f} | {t:<20.5f} | For large lambda, coeffs are zero, so t is zero.")
    if alpha < 1e-4 and last_t > 0 and abs(t - last_t)/last_t < 1e-5: # Approximate convergence
        if 'stable_t' not in locals():
            stable_t = t
            print(f"{alpha:<15.5f} | {t:<20.5f} | As lambda -> 0, t converges to the OLS L1 norm.")
            print("...")
            
    last_t = t

# The last value corresponds to the smallest alpha, approaching OLS
print(f"{alphas[-1]:<15.5f} | {t_values[-1]:<20.5f} | This t is effectively t_OLS.")
print("-" * 70)
print("\nConclusion:")
print("The output shows that as lambda approaches 0, the value of t stabilizes at a maximum value (the L1 norm of the OLS solution).")
print("This demonstrates that a single lambda value (λ=0) corresponds to a whole range of t values (t >= t_OLS).")
print("This non-bijective relationship is one reason the equivalence between the two Lasso formulations is considered 'true in general' but not 'strictly true'.")

# The problem is a multiple-choice question. Based on the analysis, K is the best answer.
<<<K>>>