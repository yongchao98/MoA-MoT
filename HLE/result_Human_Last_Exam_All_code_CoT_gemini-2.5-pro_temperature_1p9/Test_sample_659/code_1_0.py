# --- E-Score Calculation ---
e_score = 10.0
print(f"Execution (E) Score: {e_score:.1f} (assuming perfect execution)")
print("-" * 30)

# --- D-Score Calculation ---
print("Difficulty (D) Score Calculation:")

# 1. Difficulty Value (DV) from top 3 Group I elements
# Note: A 4th element is performed to avoid a short-exercise penalty but does not count for DV.
element_1 = 0.4  # D-value handstand press
element_2 = 0.4  # Another D-value handstand press
element_3 = 0.3  # C-value handstand press
total_dv = element_1 + element_2 + element_3
print(f"  Element 1 (D): {element_1:.1f}")
print(f"  Element 2 (D): {element_2:.1f}")
print(f"  Element 3 (C): {element_3:.1f}")
print(f"  Total Difficulty Value (DV): {total_dv:.1f}")

# 2. Element Group Requirements (EGR)
# Only Group I (non-acrobatic) is fulfilled. Groups II, III, and IV are not.
egr_points = 0.5
print(f"  Element Group Requirement (EGR) Points: {egr_points:.1f} (for Group I)")

# 3. Connection Value (CV)
# From connecting two D-value non-acrobatic skills (D+D >= C+C).
connection_value = 0.2
print(f"  Connection Value (CV): {connection_value:.1f} (for D+D connection)")

# Total D-Score
d_score = total_dv + egr_points + connection_value
print(f"  Total D-Score = {total_dv:.1f} (DV) + {egr_points:.1f} (EGR) + {connection_value:.1f} (CV) = {d_score:.1f}")
print("-" * 30)

# --- Neutral Deductions ---
print("Neutral Deductions:")
# Penalty for missing acrobatic element groups II and III
missing_egr_2_penalty = 0.5
missing_egr_3_penalty = 0.5
total_deductions = missing_egr_2_penalty + missing_egr_3_penalty
print(f"  Missing EGR II Penalty: {missing_egr_2_penalty:.1f}")
print(f"  Missing EGR III Penalty: {missing_egr_3_penalty:.1f}")
print(f"  Total Deductions: {total_deductions:.1f}")
print("-" * 30)

# --- Final Score Calculation ---
final_score = d_score + e_score - total_deductions
print("Final Score Calculation:")
print(f"({d_score:.1f} D-Score + {e_score:.1f} E-Score) - {total_deductions:.1f} Deductions = {final_score:.2f}")
