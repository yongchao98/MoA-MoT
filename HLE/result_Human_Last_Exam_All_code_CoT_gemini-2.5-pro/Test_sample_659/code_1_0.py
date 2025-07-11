import sys
# Redirect stdout to a variable to suppress print statements in the final output
# This is a workaround for the platform's constraints, not part of the logic.
original_stdout = sys.stdout
sys.stdout = None

# Plan:
# 1. Define the components of the gymnastics score: E-Score and D-Score.
# 2. Set the E-Score based on the "perfect execution" assumption.
# 3. Define the components of the D-Score: Difficulty Value (DV), Element Group Requirements (EGR), and Connection Value (CV).
# 4. Calculate DV:
#    - A gymnast's D-score includes up to 10 elements, but on Floor, a maximum of 5 elements from Group I (non-acrobatic) can be counted.
#    - Since the gymnast only performs handstands (Group I), only 5 skills will contribute.
#    - The 5 highest-valued handstand skills from the 2022 FIG Code of Points are one E-value (0.5) and four D-value (0.4) skills.
# 5. Calculate EGR:
#    - 0.5 points are awarded for each of the 4 Element Groups represented.
#    - The gymnast only performs skills from Group I (non-acrobatic), so only 1 group requirement is met.
# 6. Calculate CV:
#    - Connection Value on floor is only awarded for connections involving acrobatic elements.
#    - Since no acrobatic elements are performed, the CV is 0.
# 7. Sum the components to find the D-Score.
# 8. Add the D-Score and E-Score for the final score.
# 9. Print the entire calculation step-by-step, showing each number in the equations.

# Restore stdout
sys.stdout = original_stdout

# --- Calculation ---

# 1. E-Score (Execution)
# Based on "perfect execution".
e_score = 10.0

# 2. D-Score (Difficulty) Components

# 2a. Difficulty Value (DV)
# Based on the 2022 FIG Code, the 5 highest value handstand skills are:
# one E-value skill (0.5) and four D-value skills (0.4).
# Only 5 elements from Group I (non-acrobatic) can be included in the D-score calculation.
skill_values = [0.5, 0.4, 0.4, 0.4, 0.4]
difficulty_value = sum(skill_values)

# 2b. Element Group Requirements (EGR)
# There are 4 element groups on Floor. The routine only contains skills from Group I.
# The gymnast earns 0.5 points for each fulfilled group.
groups_fulfilled = 1
egr = groups_fulfilled * 0.5

# 2c. Connection Value (CV)
# On Floor, connection value requires connecting to or from an acrobatic element.
# Since only non-acrobatic handstands are performed, CV is 0.
connection_value = 0.0

# Calculate total D-Score
d_score = difficulty_value + egr + connection_value

# Calculate Final Score
final_score = d_score + e_score

# --- Output ---

print("Calculating the Maximum Score for a Handstand-Only Floor Routine")
print("-----------------------------------------------------------------")
print("The final score is composed of a Difficulty Score (D-Score) and an Execution Score (E-Score).")
print(f"Assuming perfect execution, the E-Score is {e_score:.1f}.")
print("\nThe D-Score is calculated from three components:")
print("\n1. Difficulty Value (DV):")
print("   - Only 5 non-acrobatic (Group I) skills can be counted on Floor.")
print(f"   - Top 5 handstand skill values: {skill_values[0]}, {skill_values[1]}, {skill_values[2]}, {skill_values[3]}, {skill_values[4]}")
print(f"   - Difficulty Value (DV) = {skill_values[0]} + {skill_values[1]} + {skill_values[2]} + {skill_values[3]} + {skill_values[4]} = {difficulty_value:.1f}")

print("\n2. Element Group Requirements (EGR):")
print("   - The gymnast only performs skills from Group I, meeting 1 of 4 group requirements.")
print(f"   - Element Group Requirements (EGR) = {groups_fulfilled} * 0.5 = {egr:.1f}")

print("\n3. Connection Value (CV):")
print("   - On Floor, connections must involve acrobatic elements. No CV is possible.")
print(f"   - Connection Value (CV) = {connection_value:.1f}")

print("\nTotal D-Score Calculation:")
print(f"D-Score = DV + EGR + CV")
print(f"D-Score = {difficulty_value:.1f} + {egr:.1f} + {connection_value:.1f} = {d_score:.1f}")

print("\nFinal Score Calculation:")
print(f"Final Score = D-Score + E-Score")
print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")

print(f"\n<<<12.6>>>")