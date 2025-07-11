import pandas as pd
import numpy as np
import statsmodels.api as sm

def run_simulation(n=50000, tau=10):
    """
    Simulates data to test the sign of the treatment effect coefficient.

    Args:
        n (int): Number of observations in the simulation.
        tau (float): The true, constant, positive treatment effect.
    """
    # Set seed for reproducibility
    np.random.seed(42)

    # 1. Define the Data Generating Process (DGP)
    # This DGP is for Cases 2 and 3 where D is assigned based on X.
    
    # Generate X from {-1, 0, 1}
    x_vals = [-1, 0, 1]
    X = np.random.choice(x_vals, n, p=[1/3., 1/3., 1/3.])

    # Y(0) is a non-linear function of X.
    # We add a small noise term.
    error_y0 = np.random.normal(0, 1, n)
    Y0 = 100 * X**2 + error_y0

    # The treatment effect tau is positive for everyone.
    Y1 = Y0 + tau

    # P(D=1|X) is a non-linear function of X (the "propensity score").
    # This matches the conditional random assignment assumption.
    # P(D=1|X=-1)=0.2, P(D=1|X=0)=0.8, P(D=1|X=1)=0.2
    propensity_scores = np.zeros(n)
    propensity_scores[X == -1] = 0.2
    propensity_scores[X == 0] = 0.8
    propensity_scores[X == 1] = 0.2
    
    # Generate treatment assignment D based on the propensity score
    D = np.random.binomial(1, p=propensity_scores)

    # Observed outcome Y depends on treatment status D
    Y = (1 - D) * Y0 + D * Y1

    # Create a pandas DataFrame
    df = pd.DataFrame({'Y': Y, 'D': D, 'X': X})
    df_with_const = sm.add_constant(df)

    print(f"Simulation with true positive treatment effect = {tau}\n")

    # 2. Case 2: Regress Y on D (omitting X)
    print("--- Case 2: Regress Y on D ---")
    model_case2 = sm.OLS(df['Y'], df_with_const[['const', 'D']])
    results_case2 = model_case2.fit()
    b0_c2, b1_c2 = results_case2.params
    print(f"Result: The coefficient on D is {b1_c2:.4f}, which can be negative due to OVB.")
    print("The estimated equation is:")
    print(f"Y = {b0_c2:.2f} + ({b1_c2:.2f})*D")
    print("\n" + "="*50 + "\n")


    # 3. Case 3: Regress Y on D and X
    print("--- Case 3: Regress Y on D and X ---")
    model_case3 = sm.OLS(df['Y'], df_with_const[['const', 'D', 'X']])
    results_case3 = model_case3.fit()
    b0_c3, b1_c3, b2_c3 = results_case3.params
    print(f"Result: The coefficient on D is {b1_c3:.4f}, which can be negative due to functional form misspecification.")
    print("The estimated equation is:")
    print(f"Y = {b0_c3:.2f} + ({b1_c3:.2f})*D + {b2_c3:.2f}*X")


# Run the simulation
# With tau=10, the misspecification bias in Case 3 is strong enough
# to make the coefficient on D negative.
run_simulation(tau=10)
