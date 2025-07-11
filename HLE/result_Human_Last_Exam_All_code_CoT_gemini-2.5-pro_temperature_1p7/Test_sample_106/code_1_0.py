import numpy as np
import statsmodels.api as sm

# Set a seed for reproducibility
np.random.seed(42)

# --- 0. Population and Potential Outcomes Setup ---
N = 100000  # Population size
# Generate X (pre-treatment income) for the population
X = np.random.gamma(shape=3, scale=20, size=N)

# Define the causal relationship: Y(0) is income without treatment and depends on X.
# We make the dependence strong to illustrate the potential for bias.
beta_X = 5.0
Y0 = 50 + beta_X * X + np.random.normal(0, 25, N)

# Define the treatment effect (TE) as a constant positive value for all individuals.
TE = 20.0
Y1 = Y0 + TE
# The observed income Y will be Y1 if treated (D=1), and Y0 if not (D=0).

print("This simulation will explore three cases to determine when the coefficient on a treatment D must be positive.")
print(f"The setup ensures the true treatment effect is always positive for every individual (true effect = +{TE:.2f}).")
print("-" * 80)


# --- Case 1: D is randomly assigned ---
print("Case 1: Treatment D is randomly assigned.")
print("We regress Y on a constant and D.")
D_case1 = np.random.binomial(1, 0.5, N)
Y_case1 = Y1 * D_case1 + Y0 * (1 - D_case1)

# Run the regression
X_model1 = sm.add_constant(D_case1)
model1 = sm.OLS(Y_case1, X_model1).fit()
coef_D_case1 = model1.params[1]

print(f"\nThe regression equation is: Y = {model1.params[0]:.2f} + {coef_D_case1:.2f}*D")
print(f"Result: The coefficient on D is {coef_D_case1:.2f}, which is positive and close to the true ATE ({TE:.2f}).")
print("This is because random assignment ensures the treated and control groups are comparable.")
print("The coefficient MUST be positive in this case.")
print("-" * 80)


# --- Case 2: D is conditional on X, but X is omitted ---
print("Case 2: D is assigned based on X (e.g., a program targeting low-income individuals),")
print("but we run a simple regression of Y on D, omitting X.")

# Treatment D is made more likely for individuals with lower X to create confounding.
prob_d_case23 = 1 / (1 + np.exp(-1 * (4 - 0.05 * X)))
D_case23 = np.random.binomial(1, prob_d_case23, N)
Y_case23 = Y1 * D_case23 + Y0 * (1 - D_case23)

# Run the regression omitting X
X_model2 = sm.add_constant(D_case23)
model2 = sm.OLS(Y_case23, X_model2).fit()
coef_D_case2 = model2.params[1]

print(f"\nThe regression equation is: Y = {model2.params[0]:.2f} + {coef_D_case2:.2f}*D")
print(f"Result: The coefficient on D is {coef_D_case2:.2f}, which is negative.")
print("This is due to Omitted Variable Bias. The program selects for lower-income individuals,")
print("who have lower Y outcomes regardless of treatment. This negative selection bias overwhelms the positive treatment effect.")
print("The coefficient is NOT guaranteed to be positive in this case.")
print("-" * 80)

# --- Case 3: D is conditional on X, and X is included in regression ---
print("Case 3: D is assigned based on X, and we correctly include X in the regression.")

# We use the same data from Case 2.
# Run the multiple regression including X
X_model3 = sm.add_constant(np.column_stack((D_case23, X)))
model3 = sm.OLS(Y_case23, X_model3).fit()
coef_D_case3 = model3.params[1]
coef_X_case3 = model3.params[2]

print(f"\nThe regression equation is: Y = {model3.params[0]:.2f} + {coef_D_case3:.2f}*D + {coef_X_case3:.2f}*X")
print(f"Result: The coefficient on D is {coef_D_case3:.2f}, which is positive and close to the true ATE ({TE:.2f}).")
print("By controlling for X, we remove the selection bias and correctly estimate the treatment effect.")
print("The coefficient MUST be positive in this case.")
print("-" * 80)

# --- Final Conclusion ---
print("\nConclusion:")
print("The coefficient on D is guaranteed to be positive only in Case 1 and Case 3.")