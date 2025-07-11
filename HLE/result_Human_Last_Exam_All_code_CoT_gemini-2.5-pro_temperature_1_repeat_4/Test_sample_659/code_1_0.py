# 1. Define the E-Score based on perfect execution.
e_score = 10.0

# 2. Define the Difficulty Value (DV) from the top 5 unique handstand skills.
# C-part: Press to handstand from straddle/pike (0.3)
# C-part: 1/1 turn on one arm in handstand (0.3)
# B-part: Press to handstand with straight arms/body (0.2)
# B-part: 1/1 turn in handstand (0.2)
# A-part: Press to handstand (0.1)
skill_1_c = 0.3
skill_2_c = 0.3
skill_3_b = 0.2
skill_4_b = 0.2
skill_5_a = 0.1
difficulty_value = skill_1_c + skill_2_c + skill_3_b + skill_4_b + skill_5_a

# 3. Define the Element Group Requirements (EGR).
# Only Group I (non-acrobatic) is fulfilled, earning 0.5 points.
element_group_requirements = 0.5

# 4. Define the Connection Bonus (CV).
# Two B+C connections are possible, each worth 0.2. Max CV for Group I connections is 0.4.
# Connection 1: Press (C) -> Turn (B)
# Connection 2: Press (B) -> Turn on one arm (C)
connection_1_cv = 0.2
connection_2_cv = 0.2
connection_bonus = connection_1_cv + connection_2_cv

# 5. Calculate the total D-Score.
d_score = difficulty_value + element_group_requirements + connection_bonus

# 6. Define the penalty for no dismount.
dismount_penalty = 0.5

# 7. Calculate the final score.
final_score = d_score + e_score - dismount_penalty

# 8. Print the breakdown of the calculation.
print("Calculating the Highest Possible Score:")
print("-" * 40)
print("Difficulty Value (DV) Calculation:")
print(f"Top 5 Skills = {skill_1_c} + {skill_2_c} + {skill_3_b} + {skill_4_b} + {skill_5_a} = {difficulty_value:.1f}")
print("\nD-Score Calculation:")
print(f"D-Score = {difficulty_value:.1f} (DV) + {element_group_requirements:.1f} (EGR) + {connection_bonus:.1f} (CV) = {d_score:.1f}")
print("\nFinal Score Calculation:")
print(f"Final Score = {d_score:.1f} (D-Score) + {e_score:.1f} (E-Score) - {dismount_penalty:.1f} (Penalty)")
print(f"Total = {final_score:.1f}")
print("-" * 40)
