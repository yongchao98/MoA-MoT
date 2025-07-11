# --- E-Score Calculation ---
# The E-Score starts at 10.0 for perfect execution.
base_e_score = 10.0
# A neutral deduction is applied because the routine consists only of
# non-acrobatic (Group I) elements.
composition_deduction = 0.5
final_e_score = base_e_score - composition_deduction

# --- D-Score Calculation ---
# The D-Score is the sum of Difficulty Value, Element Group points, and Connection Bonus.

# 1. Difficulty Value (DV) from the top 10 Group I skills
# The gymnast can perform up to 10 non-acrobatic elements.
# We select the highest value ones available.
skill_values = [
    0.5, 0.5,  # 2 'E' value skills
    0.4, 0.4, 0.4, 0.4, 0.4, 0.4,  # 6 'D' value skills
    0.3, 0.3,  # 2 'C' value skills
]
difficulty_value = sum(skill_values)

# 2. Element Group (EG) Points
# The gymnast only fulfills the requirement for Group I (non-acrobatic skills).
element_group_points = 0.5

# 3. Connection Bonus (CB)
# Connecting elements of C-value or higher gives a bonus. All 10 skills are C+.
# We can create three connection sequences of 3 elements each to maximize the bonus.
# Each sequence of 3 elements gives a 0.2 bonus.
connection_bonus = 0.2 + 0.2 + 0.2

# Calculate the total D-Score
d_score = difficulty_value + element_group_points + connection_bonus

# --- Final Score Calculation ---
total_score = d_score + final_e_score

# --- Print the results ---
print("Calculating the Highest Possible Score for a Handstand-Only Floor Routine")
print("-" * 70)
print(f"E-Score Calculation:")
print(f"  Perfect Execution Score: {base_e_score}")
print(f"  Composition Deduction (only non-acrobatic elements): -{composition_deduction}")
print(f"  Final E-Score = {base_e_score} - {composition_deduction} = {final_e_score:.1f}\n")

print(f"D-Score Calculation:")
print(f"  Difficulty Value (DV from 10 skills): {difficulty_value:.1f}")
print(f"  Element Group Points (EG for Group I): {element_group_points:.1f}")
print(f"  Connection Bonus (CB): {connection_bonus:.1f}")
print(f"  Final D-Score = {difficulty_value:.1f} (DV) + {element_group_points:.1f} (EG) + {connection_bonus:.1f} (CB) = {d_score:.1f}\n")

print("Final Score Breakdown:")
print(f"Total Score = D-Score + E-Score")
print(f"Total Score = ({difficulty_value:.1f} + {element_group_points:.1f} + {connection_bonus:.1f}) + ({base_e_score:.1f} - {composition_deduction:.1f})")
print(f"Total Score = {d_score:.1f} + {final_e_score:.1f} = {total_score:.1f}")

<<<14.6>>>