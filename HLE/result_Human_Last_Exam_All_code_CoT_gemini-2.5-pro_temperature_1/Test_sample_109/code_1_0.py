import numpy as np
import pandas as pd
from linearmodels.iv import IV2SLS

def perform_j_test_simulation():
    """
    This function simulates data for an instrumental variable scenario
    to demonstrate the use of the Hansen J-test for instrument validity.
    """
    # Set a seed for reproducibility
    np.random.seed(123)
    n_obs = 2000

    # 1. DEFINE THE DATA GENERATING PROCESS
    # U is an unobserved confounder that affects both the outcome Y and the invalid instrument Z2.
    U = np.random.normal(size=n_obs)

    # Z1 is a valid instrument, uncorrelated with U.
    Z1 = np.random.normal(size=n_obs)

    # Z2 is an invalid instrument because it is correlated with the confounder U.
    # This violates the exogeneity assumption.
    Z2 = 0.6 * U + np.random.normal(size=n_obs)

    # D is the binary treatment. Its assignment depends on both instruments.
    # This ensures the instruments are relevant (first stage).
    propensity = 1 / (1 + np.exp(-(0.8 * Z1 + 0.8 * Z2 - 0.5 * U)))
    D = (np.random.uniform(size=n_obs) < propensity).astype(int)

    # Y is the continuous outcome. It is affected by the treatment D and the confounder U.
    # The true treatment effect is 10.
    Y = 5 + 10 * D + 3 * U + np.random.normal(size=n_obs)

    # Create a pandas DataFrame
    data = pd.DataFrame({'Y': Y, 'D': D, 'Z1': Z1, 'Z2': Z2})
    data['const'] = 1

    # 2. PERFORM THE IV REGRESSION AND J-TEST
    # We use 2SLS with Y as the dependent variable, D as the endogenous variable,
    # and Z1 and Z2 as the instruments. The model is overidentified because we have
    # 2 instruments for 1 endogenous variable.
    model = IV2SLS(dependent=data.Y,
                   exog=data.const,
                   endog=data.D,
                   instruments=data[['Z1', 'Z2']])

    results = model.fit(cov_type='robust')

    # 3. PRINT AND INTERPRET THE RESULTS
    print("--- IV Regression Results ---")
    print(results)

    # Extract the J-test statistic and p-value
    j_statistic = results.sargan.stat
    p_value = results.sargan.pval

    print("\n--- Hansen J-Test for Overidentifying Restrictions ---")
    print(f"Null Hypothesis: All instruments (Z1 and Z2) are valid (exogenous).")
    print(f"J-Statistic: {j_statistic:.4f}")
    print(f"P-value: {p_value:.4f}")

    print("\n--- Conclusion ---")
    if p_value < 0.05:
        print("The p-value is less than 0.05, so we reject the null hypothesis.")
        print("Given that we know Z1 is a valid instrument, this test provides strong evidence")
        print("that the additional instrument, Z2, is invalid.")
    else:
        print("The p-value is greater than 0.05, so we fail to reject the null hypothesis.")
        print("The test does not find evidence to suggest that Z2 is an invalid instrument.")

# Run the simulation and print the output
perform_j_test_simulation()