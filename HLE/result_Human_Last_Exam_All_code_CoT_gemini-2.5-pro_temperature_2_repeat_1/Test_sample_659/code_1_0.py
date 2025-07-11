# --- Score Calculation based on 2022 FIG MAG Code of Points ---

# 1. E-Score (Execution)
# Assumes perfect execution, so no deductions from the starting 10.0.
e_score = 10.0

# 2. D-Score (Difficulty)
# The D-Score is composed of Difficulty Value (DV), Element Group (EG), and Connection Value (CV).

# 2a. Difficulty Value (DV)
# The top 5 non-acrobatic skills (handstands) are chosen.
# D-value skills are worth 0.4 and C-value skills are worth 0.3.
# Top 5 handstand skills: D, D, D, C, C
skill_1_value = 0.4  # Press to Handstand (D)
skill_2_value = 0.4  # Press to Handstand (D)
skill_3_value = 0.4  # Handstand 540-degree turn (D)
skill_4_value = 0.3  # Press to Handstand (C)
skill_5_value = 0.3  # Handstand 360-degree turn (C)
dv_score = skill_1_value + skill_2_value + skill_3_value + skill_4_value + skill_5_value

# 2b. Element Group (EG) Requirements
# The gymnast only performs skills from Group I (Non-Acrobatic).
# This fulfills one of the four groups, earning 0.5 points.
eg_score = 0.5

# 2c. Connection Value (CV)
# Bonus from connecting a Press (D), to a 360 Turn (C), to a 540 Turn (D).
# Connection D + C = 0.2
# Connection C + D = 0.2
cv_score = 0.2 + 0.2

# Calculate the total D-Score
d_score = dv_score + eg_score + cv_score

# Calculate the final total score
final_score = e_score + d_score

# --- Final Output ---
print("Calculating the Highest Possible Score for a Handstand-Only Routine")
print("--------------------------------------------------------------------")
print(f"E-Score (Perfect Execution): {e_score}")
print("\nD-Score is composed of three parts:")
print(f"1. Difficulty Value (Top 5 Skills): {skill_1_value} + {skill_2_value} + {skill_3_value} + {skill_4_value} + {skill_5_value} = {dv_score:.1f}")
print(f"2. Element Group Requirement (Group I only): {eg_score}")
print(f"3. Connection Value (D+C+D Series Bonus): {cv_score:.1f}")
print("--------------------------------------------------------------------")
print(f"Total D-Score Equation: {dv_score:.1f} (DV) + {eg_score:.1f} (EG) + {cv_score:.1f} (CV) = {d_score:.1f}")
print(f"Final Score Equation: {e_score:.1f} (E-Score) + {d_score:.1f} (D-Score) = {final_score:.1f}")

<<<12.7>>>