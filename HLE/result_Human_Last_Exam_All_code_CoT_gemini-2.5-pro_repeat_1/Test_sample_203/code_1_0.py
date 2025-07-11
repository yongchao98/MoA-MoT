import math

# Define the given values
altitude_AD = 12
angle_bisector_AE = 13

# --- Step 1: Calculate the lower bound for m ---
# The lower bound is derived from the fact that m > AE
lower_bound_m = angle_bisector_AE

# --- Step 2: Calculate the upper bound for m ---
# The upper bound is derived from the condition that angle A is acute.
# The derivation shows m < 156/7.
upper_bound_m_num = 156
upper_bound_m_den = 7
upper_bound_m_float = upper_bound_m_num / upper_bound_m_den

# --- Step 3: Output the final result ---
# The problem requires printing the numbers in the final equation.
# The "final equation" is the inequality defining the range of m.
print("The range of values for m for which angle A is acute is given by the inequality:")
print(f"{lower_bound_m} < m < {upper_bound_m_num}/{upper_bound_m_den}")
print(f"(Approximately, {lower_bound_m} < m < {upper_bound_m_float:.4f})")
