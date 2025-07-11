# 1. Define Execution Score (E-score)
# Assumed to be perfect as per the problem description.
e_score = 10.0

# 2. Define Difficulty Score (D-score) components
# 2.1. Difficulty Value (DV) of the highest 5 allowed Group I skills
# The routine is limited to 5 non-acrobatic (Group I) skills, with a max of 3 holds.
# We select the highest value skills that meet these criteria.
# A = 0.1, B = 0.2, C = 0.3, D = 0.4, E = 0.5

# Highest value handstand-related skills (3 holds, 2 non-holds):
skill_1 = {"name": "Planche (2s) press to handstand (2s)", "value": 0.4} # D-value, Hold
skill_2 = {"name": "V-sit (2s) press to handstand (2s)", "value": 0.3}     # C-value, Hold
skill_3 = {"name": "One-arm handstand (2s)", "value": 0.3}               # C-value, Hold
skill_4 = {"name": "Nakayama (Handstand hop 360 turn)", "value": 0.3}  # C-value, Non-Hold
skill_5 = {"name": "Guczoghy (Flair to handstand)", "value": 0.3}         # C-value, Non-Hold

skills = [skill_1, skill_2, skill_3, skill_4, skill_5]
dv_score = sum(skill['value'] for skill in skills)

# 2.2. Element Group (EG) Bonus
# The gymnast receives 0.5 points for each of the 4 Element Groups fulfilled.
# Only Group I (non-acrobatic) is fulfilled.
eg_bonus = 0.5

# 2.3. Connection Value (CV)
# CV on floor is only awarded for acrobatic connections.
cv_bonus = 0.0

# Calculate the total D-score
d_score = dv_score + eg_bonus + cv_bonus

# 3. Calculate Final Score
final_score = d_score + e_score

# --- Output the results ---
print("Calculating the Maximum Score for a Handstand-Only Floor Routine")
print("-" * 60)
print(f"E-Score (Perfect Execution): {e_score:.1f}")
print("\n--- D-Score Calculation ---")

print("\nDifficulty Value (DV) - Top 5 allowable skills from Group I:")
for i, skill in enumerate(skills):
    print(f"  Skill {i+1} ({skill['name']}): {skill['value']:.1f}")
print(f"Total Difficulty Value (DV): {skill_1['value']:.1f} + {skill_2['value']:.1f} + {skill_3['value']:.1f} + {skill_4['value']:.1f} + {skill_5['value']:.1f} = {dv_score:.1f}")

print(f"\nElement Group (EG) Bonus (Only Group I fulfilled): {eg_bonus:.1f}")
print(f"Connection Value (CV) Bonus (No acrobatic connections): {cv_bonus:.1f}")

print("\nTotal D-Score = DV + EG + CV")
print(f"Total D-Score = {dv_score:.1f} + {eg_bonus:.1f} + {cv_bonus:.1f} = {d_score:.1f}")

print("\n--- Final Score Calculation ---")
print("Final Score = D-Score + E-Score")
print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")
print("-" * 60)