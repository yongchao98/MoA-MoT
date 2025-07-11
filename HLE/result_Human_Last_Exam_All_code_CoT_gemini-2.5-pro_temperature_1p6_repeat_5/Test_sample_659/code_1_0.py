def calculate_max_handstand_score():
    """
    Calculates the highest possible score for a MAG floor routine
    consisting of only handstand skills, assuming perfect execution
    and available connections under the 2022 FIG Code of Points.
    """

    # 1. E-Score (Execution)
    # Assumed to be perfect.
    e_score = 10.0

    # 2. D-Score (Difficulty)
    # The D-Score is composed of Difficulty Value (DV), Element Group (EG) bonus,
    # and Connection Value (CV) bonus.

    # 2a. Difficulty Value (DV): Sum of the 9 best skills.
    # We select the 9 highest-rated, distinct handstand skills from the CoP.
    # D-skills (0.4 each)
    skill_1 = 0.4  # V-sit press with straight arms to handstand
    skill_2 = 0.4  # Handstand with 720 deg. turn on one arm
    skill_3 = 0.4  # Press from Japanese Swallow to handstand
    # C-skills (0.3 each)
    skill_4 = 0.3  # Press to handstand with straight arms/legs together
    skill_5 = 0.3  # V-sit press to handstand
    skill_6 = 0.3  # Handstand with 360 deg. turn on one arm
    skill_7 = 0.3  # Straddle L-sit to Stalder press to handstand
    # B-skills (0.2 each)
    skill_8 = 0.2  # Felge to handstand
    skill_9 = 0.2  # Press from straddle L-sit to handstand

    difficulty_skills = [skill_1, skill_2, skill_3, skill_4, skill_5, skill_6, skill_7, skill_8, skill_9]
    difficulty_value = sum(difficulty_skills)

    # 2b. Element Group (EG) Bonus: 0.5 for each of 4 groups fulfilled.
    # EG I (Non-acrobatic): Yes (+0.5)
    # EG II (Acrobatic forward): No (+0.0)
    # EG III (Acrobatic backward): No (+0.0)
    # EG IV (Dismount - must be acrobatic): No (+0.0)
    element_group_bonus = 0.5

    # 2c. Connection Value (CV) Bonus: Bonus for connecting skills.
    # We can maximize this by performing 4 connections with the available skills.
    # D(I) + D(I) = 0.2
    # C(I) + D(I) = 0.2
    # C(I) + C(I) = 0.1
    # B(I) + C(I) = 0.1
    connection_value_bonus = 0.2 + 0.2 + 0.1 + 0.1

    # Total D-Score
    d_score = difficulty_value + element_group_bonus + connection_value_bonus

    # 3. Penalties
    # A 0.5 penalty is applied by the Chief Judge for lack of an acrobatic dismount.
    penalty = 0.5

    # 4. Final Score Calculation
    final_score = d_score + e_score - penalty

    # --- Output the results ---
    print("Calculating the Maximum MAG Floor Score with Handstands Only\n")
    print("D-Score Breakdown:")
    print(f"  - Difficulty Value (DV) from 9 skills: {difficulty_value:.1f}")
    print(f"  - Element Group (EG) Bonus: {element_group_bonus:.1f}")
    print(f"  - Connection Value (CV) Bonus: {connection_value_bonus:.1f}")
    print(f"Total D-Score: {d_score:.1f}\n")

    print("E-Score (perfect execution): 10.0\n")
    print(f"Penalty (no acrobatic dismount): -{penalty:.1f}\n")

    # Build the final equation string with all numbers
    skill_values_str = " + ".join(map(str, difficulty_skills))
    
    print("Final Score Calculation:")
    final_equation = (f"Final Score = (DV: {skill_values_str}) + "
                      f"(EG: {element_group_bonus:.1f}) + "
                      f"(CV: {connection_value_bonus:.1f}) + "
                      f"(E-Score: {e_score:.1f}) - "
                      f"(Penalty: {penalty:.1f})")
                      
    print(final_equation)
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} - {penalty:.1f} = {final_score:.1f}")

calculate_max_handstand_score()