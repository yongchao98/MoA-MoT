def calculate_gymnast_score():
    """
    Calculates the highest possible score for a MAG floor routine
    containing only handstand skills, assuming perfect execution.
    """

    # 1. E-Score (Execution)
    # The gymnast has perfect execution.
    e_score = 10.0
    print("--- E-Score Calculation ---")
    print(f"Perfect Execution Score: {e_score}\n")

    # 2. D-Score (Difficulty)
    print("--- D-Score Calculation ---")

    # 2a. Difficulty Value (DV) from top 10 skills
    # Based on the 2022 FIG Code of Points for Group I elements.
    # We select the highest value distinct handstand skills available.
    skills = {
        "Honma press to handstand (D)": 0.4,
        "Press from L-sit to handstand (D)": 0.4,
        "Nakayama press to handstand (C)": 0.3,
        "Press from straddle L-sit to handstand (C)": 0.3,
        "Press to handstand, bent arms/straight body (B)": 0.2,
        "Press to handstand, bent arms/straddled (B)": 0.2,
        "Roll forward to handstand (A)": 0.1,
        "Handstand hold, 2 sec (A)": 0.1,
        "Japanese handstand hold, 2 sec (A)": 0.1,
        "Straddled handstand hold, 2 sec (A)": 0.1,
    }

    print("Difficulty Value (DV) is the sum of the top 10 skills:")
    dv_total = 0
    for skill, value in skills.items():
        print(f"- {skill}: {value}")
        dv_total += value
    dv_total = round(dv_total, 1)
    print(f"Total Difficulty Value (DV): {dv_total}\n")

    # 2b. Element Group Requirements (EGR)
    # 4 groups worth 0.5 each. Only Group I (non-acrobatic) is fulfilled.
    egr = 0.5
    print("Element Group Requirements (EGR):")
    print("Group I (Non-Acrobatic): 0.5")
    print("Group II (Acro Forward): 0.0")
    print("Group III (Acro Backward): 0.0")
    print("Group IV (Dismount): 0.0")
    print(f"Total EGR: {egr}\n")

    # 2c. Connection Value (CV)
    # Assumes connections are made between non-acrobatic elements.
    # D+D connection = 0.2 CV
    # C+C connection = 0.1 CV
    cv = 0.2 + 0.1
    print("Connection Value (CV):")
    print("D+D Connection: 0.2")
    print("C+C Connection: 0.1")
    print(f"Total CV: {cv}\n")

    # Total D-Score
    d_score = dv_total + egr + cv
    print(f"Total D-Score = {dv_total} (DV) + {egr} (EGR) + {cv} (CV) = {d_score}\n")

    # 3. Penalties
    # The routine lacks a required acrobatic dismount.
    penalties = 0.5
    print("--- Penalties ---")
    print(f"No Acrobatic Dismount Penalty: {penalties}\n")

    # 4. Final Score
    final_score = d_score + e_score - penalties
    print("--- Final Score Calculation ---")
    print(f"D-Score + E-Score - Penalties")
    print(f"{d_score} + {e_score} - {penalties} = {final_score}")


calculate_gymnast_score()
<<<12.5>>>