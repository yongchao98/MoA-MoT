# --- Score Calculation Breakdown ---

# 1. E-Score (Execution)
# Assumed to be perfect.
e_score = 10.0

# 2. D-Score (Difficulty)
# The D-Score is the sum of Skill Values, Element Group bonuses, and Connection Bonuses.

# 2.1. Skill Values
# A gymnast can count up to 8 acrobatic and 2 non-acrobatic skills.
# Since only handstands (non-acrobatic) are performed, only 2 skills count.
# The highest value for a handstand skill is an 'E' part, worth 0.5.
# e.g., "Press from V-sit to handstand, hold 2 sec."
skill_1_value = 0.5
skill_2_value = 0.5
total_skill_value = skill_1_value + skill_2_value

# 2.2. Element Group (EG) Requirements
# Each of the 4 groups fulfilled adds 0.5 to the D-Score.
# EG I (Non-acrobatic) is fulfilled by performing a handstand.
# EG II, III, and IV require acrobatic skills and are not fulfilled.
eg_I_bonus = 0.5
eg_II_bonus = 0.0
eg_III_bonus = 0.0
eg_IV_bonus = 0.0
total_eg_bonus = eg_I_bonus + eg_II_bonus + eg_III_bonus + eg_IV_bonus

# 2.3. Connection Bonus (CV)
# A connection of 3 Group I elements of 'D' value or higher earns 0.2.
# We assume this connection is made to maximize the score.
connection_bonus = 0.2

# Calculate the total D-Score
d_score = total_skill_value + total_eg_bonus + connection_bonus

# 3. Final Score
# Final Score = D-Score + E-Score
final_score = d_score + e_score

# --- Final Output ---
print("Calculating the maximum possible score for a floor routine with only handstand skills:")
print("\n--- D-Score Calculation ---")
print(f"Value of Top 2 Skills: {skill_1_value} + {skill_2_value} = {total_skill_value:.1f}")
print(f"Element Group Bonus (Group I only): {total_eg_bonus:.1f}")
print(f"Connection Bonus (3 connected skills): {connection_bonus:.1f}")
print(f"Total D-Score = {total_skill_value:.1f} (Skills) + {total_eg_bonus:.1f} (EG) + {connection_bonus:.1f} (CV) = {d_score:.1f}")

print("\n--- Final Score Calculation ---")
print(f"D-Score ({d_score:.1f}) + E-Score ({e_score:.1f})")
print(f"Final Score = {final_score:.1f}")
