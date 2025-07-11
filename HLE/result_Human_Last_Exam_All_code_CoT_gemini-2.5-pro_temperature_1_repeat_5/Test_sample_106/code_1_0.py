import numpy as np
import pandas as pd
import statsmodels.api as sm

def run_simulation():
    """
    Simulates a population to demonstrate treatment effect estimation
    under different regression models.
    """
    # Set up the simulation parameters
    np.random.seed(42)
    num_individuals = 100000
    
    # --- Data Generation ---
    # 1. Create the control variable X (pre-treatment income)
    # Half the population has low income (10), half has high income (20)
    pre_treatment_income = np.random.choice([10, 20], num_individuals, p=[0.5, 0.5])
    
    # 2. Assign treatment D conditional on X
    # The program is targeted at the low-income group
    prob_treatment = np.where(pre_treatment_income == 10, 0.8, 0.2)
    treatment = np.random.binomial(1, prob_treatment)
    
    # 3. Define the potential outcomes
    # The treatment effect is a +$1 increase in income for everyone.
    # This satisfies the condition that the treatment effect is always positive.
    potential_outcome_0 = pre_treatment_income  # Income without program
    potential_outcome_1 = pre_treatment_income + 1 # Income with program
    
    # 4. Observe the actual outcome Y based on treatment status
    observed_income = potential_outcome_1 * treatment + potential_outcome_0 * (1 - treatment)
    
    # Create a DataFrame for analysis
    df = pd.DataFrame({
        'Y': observed_income,
        'D': treatment,
        'X': pre_treatment_income
    })
    
    # --- Analysis ---
    print("--- Analysis of Treatment Effect Coefficients ---")

    # Case 2: D is assigned conditional on X. Regress Y on D only.
    # We expect a negative coefficient due to omitted variable bias.
    X_case2 = sm.add_constant(df['D'])
    model_case2 = sm.OLS(df['Y'], X_case2).fit()
    coef_case2 = model_case2.params['D']
    
    print("\nCase 2: Regress Y on D (omitting X)")
    print(f"The estimated coefficient on the treatment D is: {coef_case2:.4f}")
    print("This coefficient is negative because the program targets low-income individuals,")
    print("creating a strong negative selection bias that outweighs the positive treatment effect.")

    # Case 3: D is assigned conditional on X. Regress Y on D and X.
    # We expect a positive coefficient because we control for the confounder.
    X_case3 = sm.add_constant(df[['D', 'X']])
    model_case3 = sm.OLS(df['Y'], X_case3).fit()
    coef_case3 = model_case3.params['D']
    
    print("\nCase 3: Regress Y on D and X")
    print(f"The estimated coefficient on the treatment D is: {coef_case3:.4f}")
    print("This coefficient is positive, correctly identifying the direction of the treatment effect")
    print("because we have controlled for the confounding pre-treatment income X.")
    print("\nNote: In this specific data generating process (Y = X + D), the true coefficient is exactly 1.0.")


# Run the simulation and print the results
run_simulation()