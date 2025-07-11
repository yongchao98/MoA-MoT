# Plan:
# 1. Define the E-Score (Execution), which is perfect (10.0).
# 2. Define the components of the D-Score (Difficulty).
#    a. Sum of the top 9 skills: To maximize, we select four 'E' skills (0.5 value) and five 'D' skills (0.4 value).
#       - Examples of 'E' skills: Nakayama, 1080° pivot in handstand.
#       - Examples of 'D' skills: Slow roll to handstand, Guczoghy, Press from Japanese handstand.
#    b. Element Group (EG) bonus: Only EG I is fulfilled, so 0.5 points.
#    c. Connection Value (CV): Maximum is 4 connections. Connecting an E-skill to a D-skill or higher gives +0.2.
#       So, 4 * 0.2 = 0.8 points.
# 3. Sum all components to get the final score and print the equation.

# E-Score (Execution)
e_score = 10.0

# D-Score (Difficulty) components
top_9_skills = [0.5, 0.5, 0.5, 0.5, 0.4, 0.4, 0.4, 0.4, 0.4]
sum_of_skills = sum(top_9_skills)

element_group_bonus = 0.5
connection_bonus = 0.8 # 4 connections * 0.2 bonus each

# Calculate the D-Score and the Total Score
d_score = sum_of_skills + element_group_bonus + connection_bonus
total_score = d_score + e_score

# Final output showing the full equation as requested
# The final equation shows: (Sum of Skills) + EG Bonus + Connection Bonus + E-Score = Total Score
skills_str = " + ".join(map(str, top_9_skills))
print(f"The highest possible score is calculated as follows:")
print(f"Equation: (Sum of Top 9 Skills) + Element Group Bonus + Connection Bonus + E-Score")
print(f"Final Calculation: ({skills_str}) + {element_group_bonus} + {connection_bonus} + {e_score} = {total_score:.1f}")
