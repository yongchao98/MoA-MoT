import numpy as np
import statsmodels.api as sm

# --- Plan ---
# We will analyze the three cases by simulating a population and running the specified regressions.
# We will create a population where the treatment effect of D on Y is always positive for every individual.
# X will be pre-program income, D will be treatment, Y will be post-program income.
# We will check the sign of the coefficient on D in each of the three regression models.

# --- Simulation Setup ---
# 1. Define Population and Potential Outcomes
np.random.seed(42) # for reproducibility
N = 100000 # Population size

# X = pre-program income. Let's use a log-normal distribution to be realistic (income is skewed).
X = np.random.lognormal(mean=10, sigma=0.5, size=N)

# The individual treatment effect (TE) must be positive for everyone.
# We'll set it to a constant value of 2000.
TE = 2000

# Define potential outcomes Y(0) and Y(1)
# Y(0): post-program income without treatment. It depends on pre-program income (X).
u0 = np.random.normal(0, 500, size=N)
Y0 = 5000 + 1.1 * X + u0
# Y(1): post-program income with treatment.
Y1 = Y0 + TE

print("--- Simulation Setup ---")
print("This script simulates the three cases to determine when the coefficient on treatment 'D' must be positive.")
print(f"The true individual treatment effect is Y(1) - Y(0) = {TE}, which is always positive.\n")


# --- Case 1: D is randomly assigned ---
print("--- Case 1: D is randomly assigned. Regression: Y ~ constant + D ---")
# Assign treatment D completely at random (like a coin flip)
D_case1 = np.random.binomial(1, 0.5, size=N)
# Determine the observed outcome Y based on the random assignment
Y_obs_case1 = Y1 * D_case1 + Y0 * (1 - D_case1)

# Run the regression: Y on a constant and D
X_reg1 = sm.add_constant(D_case1, prepend=True)
model1 = sm.OLS(Y_obs_case1, X_reg1).fit()
b0, b1_D = model1.params

print("Population Regression Equation:")
print(f"Y = {b0:.2f} + {b1_D:.2f}*D")
print(f"The coefficient on D is {b1_D:.2f}.")
print("Conclusion: Positive. With random assignment, the coefficient on D estimates the Average Treatment Effect (ATE). Since the individual effect is always positive, the average must be positive.\n")


# --- Case 2 & 3 Setup: D assigned based on X ---
# We will now assign treatment D based on pre-program income X.
# Specifically, we'll make it more likely for individuals with lower income to receive the treatment.
# This introduces a correlation between X and D, making X a confounder.
logit_prob = 10 - 0.0003 * X
prob_d = 1 / (1 + np.exp(-logit_prob))
D_case23 = np.random.binomial(1, prob_d, size=N)
Y_obs_case23 = Y1 * D_case23 + Y0 * (1 - D_case23)


# --- Case 2: D is conditional on X, but X is omitted from regression ---
print("--- Case 2: D is conditional on X. Regression: Y ~ constant + D ---")
print("Here, we fail to control for the confounder X.")

# Run the regression: Y on a constant and D
X_reg2 = sm.add_constant(D_case23, prepend=True)
model2 = sm.OLS(Y_obs_case23, X_reg2).fit()
b0, b1_D = model2.params

print("Population Regression Equation:")
print(f"Y = {b0:.2f} + ({b1_D:.2f})*D")
print(f"The coefficient on D is {b1_D:.2f}.")
print("Conclusion: Not necessarily positive. In this simulation, it's negative due to omitted variable bias. People with low pre-program income (X) are selected for treatment, and their income (Y) would have been lower anyway. This negative selection effect can overwhelm the positive treatment effect.\n")


# --- Case 3: D is conditional on X, and X is included in regression ---
print("--- Case 3: D is conditional on X. Regression: Y ~ constant + D + X ---")
print("Here, we correctly control for the confounder X.")

# Run the regression: Y on a constant, D, and X
X_reg3 = sm.add_constant(np.column_stack((D_case23, X)), prepend=True)
model3 = sm.OLS(Y_obs_case23, X_reg3).fit()
b0, b1_D, b2_X = model3.params

print("Population Regression Equation:")
print(f"Y = {b0:.2f} + {b1_D:.2f}*D + {b2_X:.2f}*X")
print(f"The coefficient on D is {b1_D:.2f}.")
print("Conclusion: Positive. By including X in the regression, we control for the confounding variable. The coefficient on D now gives an unbiased estimate of the treatment effect, which is positive.\n")


# --- Final Answer ---
print("Based on the analysis, the coefficient on D must be positive in Case 1 and Case 3.")
