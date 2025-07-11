import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.preprocessing import StandardScaler

def solve_task():
    """
    This function replicates the R code's data generation and model comparison
    to find the expected p-value from the ANOVA test.
    """
    # Set a seed for reproducibility, similar to how R would have a default seed
    np.random.seed(123)

    # 1. Data Generation
    n = 100
    # Replicate R's recycling rule for water_availability
    water_availability = np.tile(np.array([1, 2]), n // 2)

    df = pd.DataFrame({
        'water_availability': water_availability,
        'plant_intra_div': np.random.randint(1, 3, n),
        'plant_inter_div': np.random.randint(1, 11, n)
    })

    # Simulate dependent variables based on the R script's formulas
    df['insect_richness'] = (df['water_availability'] * 0.01 +
                             df['plant_intra_div'] * 0.5 +
                             df['plant_inter_div'] * 1.2 +
                             np.random.randn(n))

    df['herbivory'] = (df['insect_richness'] * 3.14 +
                       df['water_availability'] * 0.5 +
                       df['plant_intra_div'] * 0.1 +
                       df['plant_inter_div'] * 0.2 +
                       np.random.randn(n))

    # 2. Scale the data, as done in the R script
    scaler = StandardScaler()
    df_scaled = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

    # 3. Define and fit the two models using Ordinary Least Squares (OLS)
    # Model 1: The full model
    model_1 = smf.ols('herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div', data=df_scaled).fit()

    # Model 2: The nested model (without insect_richness)
    model_2 = smf.ols('herbivory ~ water_availability + plant_intra_div + plant_inter_div', data=df_scaled).fit()

    # 4. Perform an ANOVA to compare the nested models
    anova_results = sm.stats.anova_lm(model_2, model_1)

    # The ANOVA result is a table. We will print the entire table
    # for clarity. The value we are interested in is 'Pr(>F)' in the second row.
    print("ANOVA Test Results for Model Comparison:")
    print(anova_results)

    # Extract and print the specific p-value
    p_value = anova_results.loc[1, 'Pr(>F)']
    print(f"\nThe resulting P-value from the comparison is: {p_value}")
    print("\nThis p-value is extremely close to zero, supporting the conclusion.")

solve_task()