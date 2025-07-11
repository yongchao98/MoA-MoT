import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
from scipy.stats import chi2

def solve_task():
    """
    Simulates the dataset and performs a likelihood ratio test equivalent to
    the lavaan anova() function to find the expected P-value.
    """
    # For reproducibility, same as R's default seed
    np.random.seed(42)

    n = 100
    
    # 1. Create the initial dataframe
    # In R, `1:2` used in a dataframe context with a larger n will recycle the vector.
    water_availability_raw = np.resize(np.array([1, 2]), n)
    df = pd.DataFrame({
        'water_availability': water_availability_raw,
        'plant_intra_div': np.random.randint(1, 3, size=n),
        'plant_inter_div': np.random.randint(1, 11, size=n)
    })

    # 2. Simulate the dependent variables based on the R script's logic
    # Note: R's `within` evaluates the expressions sequentially.
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra_div'] * 0.5 +
                             df['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

    # This is the key line: herbivory is strongly dependent on insect_richness
    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra_div'] * 0.1 +
                       df['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

    # 3. Scale the data as done in the R script
    scaler = StandardScaler()
    scaled_df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

    # 4. Define and fit the two nested models using Ordinary Least Squares (OLS)
    # The anova() test on these lavaan models is equivalent to a likelihood
    # ratio test on the nested OLS models for the `herbivory` outcome.

    # Model 1 (full model) for herbivory
    y = scaled_df['herbivory']
    X1 = scaled_df[['insect_richness', 'water_availability', 'plant_intra_div', 'plant_inter_div']]
    X1 = sm.add_constant(X1)
    fit_1 = sm.OLS(y, X1).fit()

    # Model 2 (restricted model) for herbivory
    X2 = scaled_df[['water_availability', 'plant_intra_div', 'plant_inter_div']]
    X2 = sm.add_constant(X2)
    fit_2 = sm.OLS(y, X2).fit()
    
    # 5. Perform the Likelihood Ratio (Chi-Square Difference) Test
    # This is what lavaan's anova() function does for nested models.
    # The test statistic is 2 * (log-likelihood of full model - log-likelihood of restricted model)
    chi_square_stat = 2 * (fit_1.llf - fit_2.llf)
    
    # The degrees of freedom is the difference in the number of parameters
    df_diff = fit_1.df_model - fit_2.df_model
    
    # The p-value is calculated from the chi-square distribution
    p_value = chi2.sf(chi_square_stat, df_diff)

    print("Likelihood Ratio Test (Chi-Square Difference Test) Results:")
    print(f"Model 1 (Full) Log-Likelihood: {fit_1.llf:.4f}")
    print(f"Model 2 (Restricted) Log-Likelihood: {fit_2.llf:.4f}")
    print("-" * 20)
    print(f"Chi-Square Difference: {chi_square_stat:.4f}")
    print(f"Degrees of Freedom Difference: {int(df_diff)}")
    # Use scientific notation for the tiny p-value
    print(f"P-value: {p_value:.2e}")
    print("\nConclusion: The p-value is extremely close to 0.")
    print("Based on the choices, the most likely reported value is 0.")

solve_task()