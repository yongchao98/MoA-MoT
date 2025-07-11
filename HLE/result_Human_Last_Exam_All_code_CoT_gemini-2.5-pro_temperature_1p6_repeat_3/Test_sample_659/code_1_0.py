#
# Calculate the maximum score for a MAG Floor routine with only handstand skills
#

# --- 1. E-Score Calculation ---
# Assuming "perfect execution" means no execution or artistry deductions.
e_score = 10.0

# --- 2. D-Score Calculation ---

# a. Element Values (EV)
# A routine's D-score is built from the top 10 elements. Handstand skills are non-acrobatic (Element Group I).
# Based on the 2022 FIG Code of Points, we select the 10 highest-value unique skills.
# We cannot get credit for both holding a position (e.g., Manna) and pressing from it.
#
# The 10 highest-value EG I skills are:
# 1. Manna press to handstand (D: 0.4)
# 2. V-sit press to handstand (C: 0.3)
# 3. One-arm handstand hold (C: 0.3)
# 4. Straddle L-sit press to handstand (B: 0.2)
# 5. Wenson, hold 2s (B: 0.2)
# 6. Press from side support to handstand (B: 0.2)
# 7. Straight arm press to handstand (A: 0.1)
# 8. Benton arm press to handstand (A: 0.1)
# 9. Handstand, hold 2s (A: 0.1)
# 10. Handstand with straddled legs, hold 2s (A: 0.1)
skill_values = [0.4, 0.3, 0.3, 0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1]
sum_of_skills = sum(skill_values)

# b. Connection Value (CV)
# On Floor, Connection Value is only awarded for acrobatic series.
# This routine has no acrobatic elements.
connection_value = 0.0

# c. Element Group Requirements (EGR)
# 0.5 bonus for each of EG I (non-acro), EG II (acro fwd), and EG III (acro bwd).
# The routine only contains EG I skills.
element_groups_fulfilled = 1
egr_bonus_per_group = 0.5
egr_bonus = element_groups_fulfilled * egr_bonus_per_group

# The total D-Score is the sum of these parts.
d_score = sum_of_skills + connection_value + egr_bonus

# --- 3. Final Score Calculation ---
final_score = d_score + e_score

# --- 4. Output the result ---
print("The final score is calculated as: D-Score + E-Score")
print("D-Score = (Sum of 10 Skills) + (Connection Value) + (Element Group Bonus)")
print("\nCalculating the final score:")
# The problem asks to output each number in the final equation.
# Final Score = (Skill Sum) + (CV) + (EGR) + (E-Score)
print(f"Final Score = {sum_of_skills:.1f} + {connection_value:.1f} + {egr_bonus:.1f} + {e_score:.1f}")
print(f"Result: {final_score:.1f}")
