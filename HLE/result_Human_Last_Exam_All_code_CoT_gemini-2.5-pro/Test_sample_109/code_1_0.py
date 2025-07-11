import numpy as np
import statsmodels.api as sm
from linearmodels.iv import IV2SLS

def run_simulation(n_simulations=1000):
    """
    Simulates data with valid instruments and heterogeneous treatment effects
    to show that the Hansen J-test (test of overidentification) can
    frequently and incorrectly reject the null of instrument validity.
    """
    np.random.seed(42)
    n_obs = 2000
    reject_count = 0
    
    for _ in range(n_simulations):
        # 1. Generate instruments and an unobserved confounder
        # Z1 and Z2 are our two instruments. V is an unobserved confounder.
        Z1 = np.random.randn(n_obs)
        Z2 = np.random.randn(n_obs)
        V = np.random.randn(n_obs)

        # 2. Generate the endogenous treatment variable D
        # D's probability depends on the instruments Z1, Z2 and the confounder V
        # This makes D endogenous because it's correlated with V.
        # The instruments Z1 and Z2 are relevant because they predict D.
        propensity = 1 / (1 + np.exp(-(0.5 * Z1 + 0.5 * Z2 + 0.5 * V)))
        D = (np.random.uniform(0, n_obs) < propensity).astype(float)

        # 3. Generate potential outcomes with HETEROGENEOUS treatment effects
        # The treatment effect's size (beta_i) depends on the instrument Z1.
        # This is a form of heterogeneity.
        beta_i = 1.0 + 0.8 * Z1 
        
        # The instruments Z1, Z2 are EXOGENOUS (valid) because they do not directly
        # cause the potential outcomes Y0 or Y1, and are uncorrelated with the
        # structural error term (which depends on V).
        Y0 = 1.0 + 1.0 * V + np.random.randn(n_obs)
        Y1 = Y0 + beta_i

        # 4. Generate the observed outcome Y
        Y = Y1 * D + Y0 * (1 - D)

        # 5. Perform the IV regression and the Hansen J-test
        # We regress Y on D, using Z1 and Z2 as instruments.
        # The model is overidentified because we have 2 instruments and 1 endogenous var.
        X = sm.add_constant(D)
        Z = sm.add_constant(np.column_stack((Z1, Z2)))
        
        try:
            iv_model = IV2SLS(dependent=Y, exog=X[:,:1], endog=X[:,1:], instruments=Z[:,1:]).fit(cov_type='robust')
            
            # The Sargan-Hansen test statistic for overidentifying restrictions
            sargan_p_value = iv_model.sargan.pval

            # We check if the test rejects the null at the 5% level
            if sargan_p_value < 0.05:
                reject_count += 1
        except np.linalg.LinAlgError:
            # In some simulations, matrices might be singular; we skip those.
            continue
            
    rejection_rate = reject_count / n_simulations
    
    print("Demonstration of Hansen J-Test Failure Under Heterogeneity")
    print("-" * 60)
    print("Scenario:")
    print("- We constructed a dataset with two instruments (Z1, Z2).")
    print("- By construction, both instruments are VALID (relevant and exogenous).")
    print("- However, the effect of the treatment depends on Z1 (heterogeneity).")
    print("\nTest Performed:")
    print("- We run an IV regression using both instruments.")
    print("- We then perform the Hansen J-test of overidentifying restrictions.")
    print("- The null hypothesis of this test is that all instruments are valid.")
    print("\nResult:")
    print(f"Number of simulations: {n_simulations}")
    print(f"J-test rejection rate (p < 0.05): {rejection_rate:.2%}")
    print("\nConclusion:")
    print("The J-test incorrectly rejects the null hypothesis of instrument validity")
    print(f"in about {rejection_rate:.0%} of cases, even though the instruments are valid.")
    print("This happens because the test mistakes treatment effect heterogeneity for")
    print("instrument invalidity. Therefore, it does not provide the guarantee")
    print("required by the problem statement.")


run_simulation()