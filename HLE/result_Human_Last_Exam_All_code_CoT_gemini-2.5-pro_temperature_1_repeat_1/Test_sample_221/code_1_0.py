import pandas as pd
import numpy as np
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler

def solve():
    """
    This function simulates the data and performs the model comparison as described
    in the R script to find the expected P-value.
    """
    # Set a random seed for reproducibility
    np.random.seed(42)
    n = 100

    # 1. Create the dataframe
    # We tile [1, 2] to get 50 of each for water_availability, similar to the R code's 1:2
    # which R's data.frame function would recycle.
    df = pd.DataFrame({
        'water_availability': np.tile([1, 2], n // 2),
        'plant_intra_div': np.random.randint(1, 3, size=n),
        'plant_inter_div': np.random.randint(1, 11, size=n)
    })

    # 2. Generate the dependent variables based on the formulas
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra_div'] * 0.5 +
                             df['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra_div'] * 0.1 +
                       df['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

    # 3. Scale the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df)
    scaled_df = pd.DataFrame(scaled_data, columns=df.columns)

    # 4. Define and fit the two models for the 'herbivory' outcome
    # This corresponds to the part of the SEM being tested by the anova.
    # Model 1 (full model)
    X1 = sm.add_constant(scaled_df[['insect_richness', 'water_availability', 'plant_intra_div', 'plant_inter_div']])
    y = scaled_df['herbivory']
    fit_1 = sm.OLS(y, X1).fit()

    # Model 2 (restricted model, omits insect_richness)
    X2 = sm.add_constant(scaled_df[['water_availability', 'plant_intra_div', 'plant_inter_div']])
    fit_2 = sm.OLS(y, X2).fit()

    # 5. Perform a likelihood ratio test to compare the models
    # This is equivalent to lavaan's anova(fit_1, fit_2) for this comparison.
    lr_statistic = 2 * (fit_1.llf - fit_2.llf)
    df_diff = fit_1.df_model - fit_2.df_model
    from scipy.stats import chi2
    p_value = chi2.sf(lr_statistic, df_diff)

    # The anova_lm function in statsmodels provides the same test as an F-test.
    # anova_results = sm.stats.anova_lm(fit_2, fit_1)
    # p_value_f_test = anova_results.iloc[1]['Pr(>F)']
    
    # The p-value will be extremely small. Python may display it in scientific notation.
    # For the purpose of answering the multiple-choice question, we are looking for a
    # value that is effectively zero.
    
    print(f"The Likelihood Ratio test compares the full model against the model without the 'insect_richness -> herbivory' path.")
    print(f"Model 1 Log-Likelihood: {fit_1.llf:.4f}")
    print(f"Model 2 Log-Likelihood: {fit_2.llf:.4f}")
    print(f"Chi-Square difference: {lr_statistic:.4f}")
    print(f"Degrees of freedom difference: {df_diff}")
    print(f"Resulting P-value: {p_value}")
    print("\nSince the p-value is a very small number (effectively zero), this indicates a highly significant difference between the models.")
    print("This means the path from 'insect_richness' to 'herbivory' is critical for explaining the data, as expected from the data generation process.")
    print("Therefore, the most appropriate answer is 0.")

solve()