def calculate_max_handstand_score():
    """
    Calculates the maximum possible score on MAG floor with only handstand skills
    under the 2022 FIG Code of Points, assuming perfect execution.
    """

    # 1. D-Score: Value of the top 10 skills
    # Based on the 2022 FIG MAG Code of Points, we select the 10 highest-value
    # unique handstand skills from Element Group I.
    # (Difficulty: B=0.2, C=0.3, D=0.4, E=0.5)
    top_10_skills = {
        "Manna press to handstand": 0.5,                  # E
        "Handstand with 3/1 turn (1080°)": 0.5,            # E
        "Handstand with 2.5/1 turn (900°)": 0.5,           # E
        "Straddle L-sit press to handstand": 0.4,          # D
        "Handstand with 2/1 turn (720°)": 0.4,             # D
        "Straight arm press with V-sit hold (2s)": 0.3,    # C
        "Handstand with 1.5/1 turn (540°)": 0.3,           # C
        "V-sit, roll backward to handstand (2s)": 0.3,     # C
        "Straight arm press to handstand": 0.2,            # B
        "Handstand with 1/1 turn (360°)": 0.2              # B
    }

    skill_values = list(top_10_skills.values())
    sum_of_skills = sum(skill_values)

    # 2. D-Score: Element Group Requirement (EGR) Bonus
    # All skills are from Group I (non-acrobatic). No skills from Groups II, III, or IV.
    # Therefore, only 1 group is fulfilled, granting a 0.5 bonus.
    egr_bonus = 0.5

    # 3. D-Score: Connection Value (CV)
    # The MAG CoP does not award Connection Value for non-acrobatic elements on floor.
    connection_value = 0.0

    # Total D-Score is the sum of the above components
    d_score = sum_of_skills + egr_bonus + connection_value

    # 4. E-Score (Execution Score)
    # The prompt assumes perfect execution, which starts at 10.0.
    # We are ignoring compositional deductions (e.g. for lack of acrobatic elements)
    # to adhere to the "perfect execution" spirit of the prompt.
    e_score = 10.0

    # 5. Final Score Calculation
    total_score = d_score + e_score

    # --- Outputting the results step-by-step ---
    print("Calculating the Maximum MAG Floor Score with Only Handstands")
    print("-" * 60)
    print("Component 1: Difficulty from Top 10 Skills")
    for skill, value in top_10_skills.items():
        print(f"- {skill}: {value:.1f}")
    
    skill_values_str = " + ".join(map(str, skill_values))
    print(f"\nSum of Skill Values = {skill_values_str} = {sum_of_skills:.1f}")
    print("-" * 60)

    print(f"Component 2: Element Group Requirement (EGR) Bonus")
    print(f"Only Group I fulfilled -> Bonus = {egr_bonus:.1f}")
    print("-" * 60)
    
    print(f"Component 3: Connection Value (CV)")
    print(f"No CV for non-acrobatic connections on floor -> Value = {connection_value:.1f}")
    print("-" * 60)

    print(f"Total D-Score = {sum_of_skills:.1f} (Skills) + {egr_bonus:.1f} (EGR) + {connection_value:.1f} (CV) = {d_score:.1f}")
    print("-" * 60)

    print(f"Total E-Score (assumed perfect execution) = {e_score:.1f}")
    print("-" * 60)

    print("\nFinal Score Calculation:")
    print(f"Total Score = D-Score + E-Score")
    print(f"Total Score = {d_score:.1f} + {e_score:.1f} = {total_score:.1f}")
    
    print("\nIn one line, the full calculation is:")
    # The final print statement shows every number involved in the sum
    print(f"Final Score = {skill_values_str} + {egr_bonus} + {connection_value} + {e_score} = {total_score:.1f}")


calculate_max_handstand_score()