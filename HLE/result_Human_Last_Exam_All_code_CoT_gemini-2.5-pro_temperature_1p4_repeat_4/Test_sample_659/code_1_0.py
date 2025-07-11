def calculate_max_handstand_score():
    """
    Calculates the maximum possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution.
    """
    # 1. E-Score (Execution)
    # The problem assumes perfect execution.
    e_score = 10.0

    # 2. D-Score (Difficulty)
    # The D-Score is the sum of Skill Values, Element Group Requirements (EGR),
    # and Connection Bonus (CV).

    # 2a. Skill Values
    # To maximize the score, we select the highest value skills and those that
    # allow for the maximum connection bonus. We need 8 skills of 'D' value
    # or higher to get the max CV of 0.8.
    
    skills = {
        "Handstand lower to Manna (F)": 0.6,
        "Manna (E)": 0.5,
        "Straddle planch (D)": 0.4,
        "V-sit press to handstand (D)": 0.4,
        "Handstand lower to straddle planch (D)": 0.4,
        "Jonasson (D)": 0.4,
        "Wendeswing fwd to handstand (D)": 0.4,
        "Straddle L-sit press to handstand (D)": 0.4,
        "Straight body press to handstand (C)": 0.3,
        "From L-sit through handstand with 1/1 turn (C)": 0.3 # Dismount
    }
    skill_value_sum = sum(skills.values())

    # 2b. Element Group Requirements (EGR)
    # Handstand skills are all Group I (non-acrobatic).
    # Only 1 of 4 groups is fulfilled.
    egr_points = 0.5

    # 2c. Connection Bonus (CV)
    # Connecting two Group I elements of D-value or higher gives 0.2 CV.
    # The routine has 8 D+ skills, allowing for 4 such connections.
    num_connections = 4
    bonus_per_connection = 0.2
    connection_bonus = num_connections * bonus_per_connection

    # Calculate Total D-Score and Final Score
    d_score = skill_value_sum + egr_points + connection_bonus
    final_score = d_score + e_score

    # Print the breakdown
    print("--- Maximum MAG Floor Score with Handstands Only ---")
    print("\nE-Score (Perfect Execution): 10.0")
    print("\nD-Score Calculation:")
    print(f"  - Top 10 Skill Values Sum: {skill_value_sum:.1f}")
    print(f"  - Element Group Requirements (Group I): {egr_points:.1f}")
    print(f"  - Connection Bonus (4 x D+ connections): {connection_bonus:.1f}")
    
    print("\nFinal Score Calculation:")
    print(f"Final Score = D-Score + E-Score")
    print(f"Final Score = (Skill Sum + EGR + CV) + E-Score")
    # The final equation with all numbers
    print(f"Final Score = ({skill_value_sum:.1f} + {egr_points:.1f} + {connection_bonus:.1f}) + {e_score:.1f}")
    print(f"Final Score = {d_score:.1f} + {e_score:.1f} = {final_score:.1f}")

calculate_max_handstand_score()
<<<15.4>>>