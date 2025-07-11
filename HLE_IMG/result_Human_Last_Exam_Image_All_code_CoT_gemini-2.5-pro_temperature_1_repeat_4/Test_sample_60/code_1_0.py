import pandas as pd
import numpy as np

# Step 1: Create the dataset
data = {
    'y':  [0.27, 0.33, 0.52, 0.79, 0.89, 0.18, 0.39, 0.61, 0.31, 0.59, 0.88, 0.90, 0.54, 0.87, 0.67],
    'x1': [0.07, 0.11, 0.28, 0.62, 0.79, 0.03, 0.15, 0.37, 0.10, 0.35, 0.78, 0.81, 0.29, 0.75, 0.45],
    'x2': [0.81, 1.00, 1.57, 2.36, 2.67, 0.55, 1.16, 1.82, 0.93, 1.76, 2.65, 2.71, 1.62, 2.60, 2.01],
    'x3': [0.22, 0.27, 0.62, 0.23, 0.33, 0.21, 0.02, 0.93, 0.86, 0.29, 0.75, 0.39, 0.53, 0.21, 0.78],
    'x4': [0.98, 1.00, 0.79, 0.03, 0.56, 0.38, 0.40, 0.37, 0.23, 0.10, 0.35, 0.26, 0.73, 0.40, 0.63],
    'x5': [0.55, 0.77, 0.99, 0.48, 0.66, 0.80, 0.09, 0.54, 0.43, 0.10, 0.59, 0.51, 0.95, 0.39, 0.65],
    'x6': [0.21, 0.25, 0.73, 0.59, 0.93, 0.87, 0.76, 0.96, 0.72, 0.10, 0.13, 0.71, 0.59, 0.75, 0.33]
}
df = pd.DataFrame(data)

# Step 2: Calculate and print the correlation matrix
correlation_matrix = df.corr()
print("Correlation Matrix:")
print(correlation_matrix)

# Step 3: Analyze the correlations to select the best set of variables
print("\nAnalysis:")
# Identify endogenous and potential exogenous variables/instruments
print("x1 and x2 are endogenous as they are functions of y.")
print(f"Correlation(x1, y) = {correlation_matrix.loc['x1', 'y']:.3f}")
print(f"Correlation(x2, y) = {correlation_matrix.loc['x2', 'y']:.3f}\n")

print("Evaluating x3, x4, x5, x6 as potential instruments or exogenous regressors.")
# Find the most exogenous variable (least correlated with y)
exog_corr = correlation_matrix['y'][['x3', 'x4', 'x5', 'x6']].abs()
best_exog_var = exog_corr.idxmin()
print(f"Variable least correlated with y (best candidate for an exogenous regressor): {best_exog_var}")
print(f"Correlation({best_exog_var}, y) = {correlation_matrix.loc[best_exog_var, 'y']:.3f}\n")

# Find the most relevant instrument for x1
# We exclude best_exog_var from the candidates for the external instrument
instrument_candidates = ['x3', 'x4', 'x5', 'x6']
relevance_corr = correlation_matrix['x1'][instrument_candidates].abs()
best_instrument = relevance_corr.idxmax()
print(f"Variable most correlated with x1 (best candidate for an external instrument): {best_instrument}")
print(f"Correlation({best_instrument}, x1) = {correlation_matrix.loc[best_instrument, 'x1']:.3f}\n")

# Check if x5 is a relevant instrument
print(f"Checking relevance of x5 for x1: Corr(x5, x1) = {correlation_matrix.loc['x5', 'x1']:.3f}. This is near zero, so x5 is not a relevant external instrument.")

# Conclusion
print("\nConclusion:")
print("A robust GMM setup would consist of:")
print(f"- The endogenous variable: x1")
print(f"- An included exogenous variable (least correlated with y): {best_exog_var}")
print(f"- An external instrument for x1 (most correlated with x1): {best_instrument}")
print(f"This gives the set of variables: {{x1, {best_instrument}, {best_exog_var}}}")
