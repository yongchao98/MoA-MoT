import pandas as pd
import numpy as np
from sklearn.preprocessing import scale
from semopy import Model
from semopy.stats import compare
import warnings

# Suppress potential convergence warnings from semopy for this example
warnings.filterwarnings("ignore", category=UserWarning)

def solve_sem_anova():
    """
    This function simulates the R code's data generation and model comparison
    to determine the expected p-value from the anova test.
    """
    # Set a seed for reproducibility. R's sample() and rnorm() will produce
    # different random numbers than numpy, but the statistical properties
    # of the resulting data will be the same, leading to the same conclusion.
    np.random.seed(123)
    n = 100

    # 1. Generate data following the R code's logic
    # In R, 1:2 is recycled to length 100. np.tile achieves the same.
    water_availability = np.tile([1, 2], n // 2)
    # R's sample(1:2, ...) is equivalent to randint(1, 3, ...) in numpy
    plant_intra_div = np.random.randint(1, 3, size=n)
    plant_inter_div = np.random.randint(1, 11, size=n)

    temp_data = {
        'water_availability': water_availability,
        'plant_intra.div': plant_intra_div,
        'plant_inter.div': plant_inter_div
    }

    insect_richness = (temp_data['water_availability'] * 0.01 +
                       temp_data['plant_intra.div'] * 0.5 +
                       temp_data['plant_inter.div'] * 1.2 +
                       np.random.randn(n))

    herbivory = (insect_richness * 3.14 +
                 temp_data['water_availability'] * 0.5 +
                 temp_data['plant_intra.div'] * 0.1 +
                 temp_data['plant_inter.div'] * 0.2 +
                 np.random.randn(n))

    temp_data['insect_richness'] = insect_richness
    temp_data['herbivory'] = herbivory

    example_unscaled = pd.DataFrame(temp_data)

    # 2. Scale the data, same as the `scale()` function in R
    example_scaled_array = scale(example_unscaled)
    example = pd.DataFrame(example_scaled_array, columns=example_unscaled.columns)
    
    # semopy can't handle dots in variable names, replace with underscores
    example.columns = [c.replace('.', '_') for c in example.columns]
    
    # 3. Define the two models
    model_1_desc = """
    herbivory ~ insect_richness + water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """
    
    model_2_desc = """
    # Path from insect_richness to herbivory is removed
    herbivory ~ water_availability + plant_intra_div + plant_inter_div
    insect_richness ~ water_availability + plant_intra_div + plant_inter_div
    """

    # 4. Fit the models
    model_1 = Model(model_1_desc)
    fit_1 = model_1.fit(example)

    model_2 = Model(model_2_desc)
    fit_2 = model_2.fit(example)
    
    # 5. Compare the models using a likelihood-ratio test (equivalent to R's anova)
    # The 'compare' function returns a pandas DataFrame with the test statistics.
    comparison_results = compare(fit_1, fit_2)
    
    # The result we care about is the p-value for the Chi-squared difference test
    p_value = comparison_results['p-value']['m1 m2']

    # The effect size (3.14) is so large that the p-value will be extremely small,
    # effectively zero. We print it out to show this.
    print(f"The Likelihood Ratio Test compares Model 1 and Model 2.")
    print(f"The null hypothesis is that the simpler model (Model 2) is sufficient.")
    print(f"In Model 2, the path from 'insect_richness' to 'herbivory' is constrained to 0.")
    print(f"However, the data was generated with a strong coefficient of 3.14 for this path.")
    print(f"Therefore, we expect to strongly reject the null hypothesis.")
    print(f"The calculated P-value from the comparison is: {p_value}")
    print(f"This value is extremely small, and for the purpose of the multiple choice question, is effectively 0.")

solve_sem_anova()