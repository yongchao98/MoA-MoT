# A script to calculate the highest possible MAG Floor score with only handstand skills.

# --- D-Score Calculation ---

# 1. Element Values: A maximum of 3 Group I (non-acrobatic) skills are counted.
# We select the three highest-value skills from the Code of Points.
element_e_value = 0.5  # Manna to handstand press, hold (2s)
element_d1_value = 0.4 # Straddle support to handstand press, hold (2s)
element_d2_value = 0.4 # Press from one arm to handstand, hold (2s)
total_element_value = element_e_value + element_d1_value + element_d2_value

# 2. Element Group (EG) Bonus: 0.5 for fulfilling Group I.
eg_bonus = 0.5

# 3. Connection Value (CV) Bonus: For a D+E and E+D series of strength elements.
cv_bonus_de = 0.2
cv_bonus_ed = 0.2
total_cv_bonus = cv_bonus_de + cv_bonus_ed

# 4. Dismount Penalty: 0.5 deduction for no acrobatic dismount.
dismount_penalty = 0.5

# Calculate the total D-Score.
d_score = total_element_value + eg_bonus + total_cv_bonus - dismount_penalty

# --- E-Score Calculation ---

# Per the prompt, execution is perfect, resulting in a full E-Score.
e_score = 10.0

# --- Final Score Calculation ---

# The final score is the sum of the D-Score and E-Score.
final_score = d_score + e_score

print("D-Score Calculation:")
print(f"({element_e_value} + {element_d1_value} + {element_d2_value}) [Elements] + {eg_bonus} [EG Bonus] + {total_cv_bonus} [CV Bonus] - {dismount_penalty} [Penalty] = {d_score}")
print("\nE-Score Calculation:")
print(f"Perfect Execution Score = {e_score}")
print("\nFinal Score (D-Score + E-Score):")
print(f"{d_score} + {e_score} = {final_score}")
