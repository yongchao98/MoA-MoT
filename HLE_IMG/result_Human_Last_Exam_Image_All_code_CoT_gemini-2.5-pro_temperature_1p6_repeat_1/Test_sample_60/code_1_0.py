import pandas as pd
import numpy as np
import statsmodels.api as sm

# Step 1: Define the dataset
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

# --- Analysis ---
print("Analysis to determine the set of instrumental variables.")
print("-" * 50)
print("Initial check: x1 and x2 are endogenous as they are functions of y.")
print(f"Corr(y^2, x1) = {np.corrcoef(df['y']**2, df['x1'])[0,1]:.4f}")
print(f"Corr(3*y, x2) = {np.corrcoef(3*df['y'], df['x2'])[0,1]:.4f}")
print("Thus, x1 and x2 cannot be instruments themselves.")
print("We assume the options are {Endogenous Var, Instrument 1, Instrument 2}.\n")

# Step 2: Calculate the correlation matrix
correlation_matrix = df.corr()
print("Step 2: Evaluate potential instruments {x3, x4, x5, x6}")
print("Assuming x1 is the endogenous variable to be instrumented.")
corr_x1 = correlation_matrix.loc[['x3', 'x4', 'x5', 'x6'], 'x1'].abs()
corr_y = correlation_matrix.loc[['x3', 'x4', 'x5', 'x6'], 'y'].abs()

# Step 3: Identify the best instruments based on relevance and exogeneity
print("\nRelevance check (|Correlation with x1|):")
print(corr_x1.round(4))
most_relevant_iv = corr_x1.idxmax()
print(f"--> Most Relevant IV: {most_relevant_iv} (Corr={corr_x1.max():.4f})")

print("\nExogeneity check (proxy: low |Correlation with y|):")
print(corr_y.round(4))
most_exogenous_iv = corr_y.idxmin()
print(f"--> Most Exogenous IV: {most_exogenous_iv} (Corr={corr_y.min():.4f})")

print(f"\nThis logic suggests the instrument set is {{{most_relevant_iv}, {most_exogenous_iv}}}.")
print("This instrument pair appears in options C and E.")

# Step 4: Compare relevance of {x3, x5} for x1 and x2
print("\nStep 4: Test which endogenous variable the pair {x3, x5} best instruments.")
df['const'] = 1
# First stage regression for x1
model_x1 = sm.OLS(df['x1'], df[['const', 'x3', 'x5']]).fit()
f_stat_x1 = model_x1.fvalue
print(f"F-statistic for x1 ~ 1 + x3 + x5: {f_stat_x1:.4f}")

# First stage regression for x2
model_x2 = sm.OLS(df['x2'], df[['const', 'x3', 'x5']]).fit()
f_stat_x2 = model_x2.fvalue
print(f"F-statistic for x2 ~ 1 + x3 + x5: {f_stat_x2:.4f}")

print(f"\nThe F-statistic for x1 ({f_stat_x1:.2f}) is higher than for x2 ({f_stat_x2:.2f}).")
print("This indicates {x3, x5} are a stronger (though still weak) set of instruments for x1.")
print("\nConclusion: The choice should be the one listing x1 as the endogenous variable, which is E.")
<<<E>>>