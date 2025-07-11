import pandas as pd
import numpy as np

# 1. Define the dataset
data = {
    'y': [0.27, 0.33, 0.52, 0.79, 0.89, 0.18, 0.39, 0.61, 0.31, 0.59, 0.88, 0.90, 0.54, 0.87, 0.67],
    'x1': [0.07, 0.11, 0.28, 0.62, 0.79, 0.03, 0.15, 0.37, 0.10, 0.35, 0.78, 0.81, 0.29, 0.75, 0.45],
    'x2': [0.81, 1.00, 1.57, 2.36, 2.67, 0.55, 1.16, 1.82, 0.93, 1.76, 2.65, 2.71, 1.62, 2.60, 2.01],
    'x3': [0.22, 0.27, 0.62, 0.23, 0.33, 0.21, 0.02, 0.93, 0.86, 0.29, 0.75, 0.39, 0.53, 0.21, 0.78],
    'x4': [0.98, 1.00, 0.79, 0.03, 0.56, 0.38, 0.40, 0.37, 0.23, 0.10, 0.35, 0.26, 0.73, 0.40, 0.63],
    'x5': [0.55, 0.77, 0.99, 0.48, 0.66, 0.80, 0.09, 0.54, 0.43, 0.10, 0.59, 0.51, 0.95, 0.39, 0.65],
    'x6': [0.21, 0.25, 0.73, 0.59, 0.93, 0.87, 0.76, 0.96, 0.72, 0.10, 0.13, 0.71, 0.59, 0.75, 0.33]
}
df = pd.DataFrame(data)

# 2. Explain the logic
print("To find the most suitable instrumental variable, we need one that is:")
print("1. Relevant: Correlated with the endogenous variable.")
print("2. Exogenous: Not correlated with the error term in the causal model.\n")

print("Analysis of variables:")
print("- x1 and x2 are functions of y (x1 ≈ y^2, x2 ≈ 3y). This makes them endogenous.")
print("- Therefore, x1 and x2 cannot be valid instruments as they violate the exogeneity condition.\n")

print("Evaluating potential instruments (x3, x4, x5, x6):")
print("We assume x3, x4, x5, and x6 are exogenous as they are described as 'random variables'.")
print("The best instrument will be the most relevant, i.e., the one with the highest correlation with the endogenous variable (represented by y).\n")

# 3. Calculate correlations with y
potential_instruments = ['x3', 'x4', 'x5', 'x6']
correlations = df[potential_instruments + ['y']].corr()['y']

print("Correlations of potential instruments with y:")
best_instrument = ''
max_corr = -1

for var in potential_instruments:
    corr_val = correlations[var]
    print(f"Correlation(y, {var}) = {corr_val:.4f}")
    if abs(corr_val) > max_corr:
        max_corr = abs(corr_val)
        best_instrument = var

# 4. State the conclusion
print(f"\nThe variable with the highest absolute correlation with y is {best_instrument}, with a correlation of {correlations[best_instrument]:.4f}.")
print(f"Therefore, {best_instrument} is the most suitable and relevant instrument for identifying causal effects on y.")
