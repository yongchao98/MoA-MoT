import numpy as np

# Step 1: Define the dataset
y = np.array([0.27, 0.33, 0.52, 0.79, 0.89, 0.18, 0.39, 0.61, 0.31, 0.59, 0.88, 0.90, 0.54, 0.87, 0.67])
x1 = np.array([0.07, 0.11, 0.28, 0.62, 0.79, 0.03, 0.15, 0.37, 0.10, 0.35, 0.78, 0.81, 0.29, 0.75, 0.45])
x3 = np.array([0.22, 0.27, 0.62, 0.23, 0.33, 0.21, 0.02, 0.93, 0.86, 0.29, 0.75, 0.39, 0.53, 0.21, 0.78])
x4 = np.array([0.98, 1.00, 0.79, 0.03, 0.56, 0.38, 0.40, 0.37, 0.23, 0.10, 0.35, 0.26, 0.73, 0.40, 0.63])
x5 = np.array([0.55, 0.77, 0.99, 0.48, 0.66, 0.80, 0.09, 0.54, 0.43, 0.10, 0.59, 0.51, 0.95, 0.39, 0.65])
x6 = np.array([0.21, 0.25, 0.73, 0.59, 0.93, 0.87, 0.76, 0.96, 0.72, 0.10, 0.13, 0.71, 0.59, 0.75, 0.33])

# Step 2: Explain the logic
print("To find the most suitable instrumental variable, we need to check the 'relevance' condition.")
print("This means finding which potential instrument (x3, x4, x5, x6) has the strongest correlation with the endogenous variable (e.g., x1).\n")
print("Variables x1 and x2 are endogenous because they are functions of y (x1=y^2, x2=3y) and cannot be instruments.\n")

# Step 3: Calculate the correlation of each potential instrument with the endogenous variable x1
corr_x3_x1 = np.corrcoef(x3, x1)[0, 1]
corr_x4_x1 = np.corrcoef(x4, x1)[0, 1]
corr_x5_x1 = np.corrcoef(x5, x1)[0, 1]
corr_x6_x1 = np.corrcoef(x6, x1)[0, 1]

# Store correlations for comparison
correlations = {
    'x3': corr_x3_x1,
    'x4': corr_x4_x1,
    'x5': corr_x5_x1,
    'x6': corr_x6_x1
}

print("Correlation with endogenous variable x1:")
for var, corr in correlations.items():
    print(f"Correlation between {var} and x1: {corr:.4f}")

# Step 4: Find the instrument with the highest absolute correlation
best_instrument = max(correlations, key=lambda k: abs(correlations[k]))
highest_corr = correlations[best_instrument]

print("\n--- Conclusion ---")
print(f"The variable with the highest absolute correlation with x1 is '{best_instrument}'.")
print(f"Its correlation is {highest_corr:.4f}, with an absolute value of {abs(highest_corr):.4f}.")
print(f"Therefore, '{best_instrument}' is the most suitable instrumental variable for identifying causal effects on y.")
