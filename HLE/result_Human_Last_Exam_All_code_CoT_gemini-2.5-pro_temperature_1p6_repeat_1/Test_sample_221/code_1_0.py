import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2

def solve_case():
    """
    This function simulates the dataset and statistical test described in the problem
    to determine the expected p-value from the anova comparison.
    """
    # Set a seed for reproducibility, same as in the R code analysis
    np.random.seed(42)

    # 1. Data Generation
    n = 100
    # In R, data.frame(v=1:2) with n=100 recycles the vector to [1, 2, 1, 2, ...]
    water_availability = np.tile([1, 2], n // 2)

    df = pd.DataFrame({
        'water_availability': water_availability,
        # In R, sample(1:2) includes 2. In Python, randint is exclusive of the high end.
        'plant_intra_div': np.random.randint(1, 3, size=n),
        'plant_inter_div': np.random.randint(1, 11, size=n)
    })

    # Generate the dependent variables based on the defined relationships
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                                   df['plant_intra_div'] * 0.5 +
                                   df['plant_inter_div'] * 1.2 +
                                   np.random.randn(n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                             df['water_availability'] * 0.5 +
                             df['plant_intra_div'] * 0.1 +
                             df['plant_inter_div'] * 0.2 +
                             np.random.randn(n))

    # 2. Scaling
    # The anova result is invariant to scaling, but we scale for faithfulness to the original script.
    # R's scale() uses a divisor of N-1 for standard deviation, so we use ddof=1.
    df_scaled = (df - df.mean()) / df.std(ddof=1)

    # 3. Model Fitting
    # The lavaan anova test compares the overall model fits. The key difference lies in the
    # herbivory equation, so we can analyze that with OLS to get the log-likelihoods.
    
    # Model 1 (Full model): herbivory ~ insect_richness + ...
    formula_1 = 'herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div'
    fit_1 = smf.ols(formula_1, data=df_scaled).fit()

    # Model 2 (Restricted model): herbivory ~ ... (without insect_richness)
    formula_2 = 'herbivory ~ water_availability + plant_intra_div + plant_inter_div'
    fit_2 = smf.ols(formula_2, data=df_scaled).fit()

    # 4. Perform Likelihood Ratio Test
    log_likelihood_1 = fit_1.llf
    log_likelihood_2 = fit_2.llf
    
    # The difference in degrees of freedom is the number of constrained parameters (1).
    df_diff = fit_1.df_model - fit_2.df_model
    
    # The Likelihood Ratio (Chi-Squared Difference) statistic calculation:
    lr_statistic = 2 * (log_likelihood_1 - log_likelihood_2)
    
    # Calculate the p-value from the chi-squared distribution
    p_value = chi2.sf(lr_statistic, df_diff)

    print("--- Likelihood Ratio Test ---")
    print(f"The final test equation is: P(Chi_Square(df={int(df_diff)}) > {lr_statistic:.4f})")
    print("\nComponent values:")
    print(f"Log-Likelihood of Full Model (fit_1): {log_likelihood_1:.4f}")
    print(f"Log-Likelihood of Restricted Model (fit_2): {log_likelihood_2:.4f}")
    print(f"Likelihood Ratio Statistic (Chi-Square value): {lr_statistic:.4f}")
    print(f"Degrees of Freedom Difference: {int(df_diff)}")
    print(f"\nFinal P-value: {p_value}")


solve_case()