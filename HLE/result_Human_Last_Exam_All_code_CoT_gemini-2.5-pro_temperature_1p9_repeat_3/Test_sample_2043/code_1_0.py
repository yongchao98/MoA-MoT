def solve_ems_scenario():
    """
    Calculates the best hospital destination based on a clinical scoring model.
    """
    # Patient condition weights
    trauma_urgency = 10
    toxicology_urgency = 5
    time_penalty_factor = 3

    # Hospital options: [Name, Trauma Capability, Toxin Capability, Transport Time (min), Letter]
    hospitals = [
        ("Level 4 trauma center", 2, 1, 6, "A"),
        ("Level 3 trauma center", 5, 1, 7, "B"),
        ("Level 2 trauma center", 8, 1, 8, "C"),
        ("Level 2 trauma center with a toxicologist", 8, 10, 15, "D"),
        ("Level 1 trauma center with toxicologist", 10, 10, 15, "E")
    ]

    best_hospital = None
    max_score = -float('inf')

    print("Analyzing hospital options based on patient need...\n")

    for name, trauma_cap, tox_cap, time, letter in hospitals:
        score = (trauma_cap * trauma_urgency) + (tox_cap * toxicology_urgency) - (time * time_penalty_factor)
        if score > max_score:
            max_score = score
            best_hospital = {
                "name": name,
                "trauma_cap": trauma_cap,
                "tox_cap": tox_cap,
                "time": time,
                "letter": letter,
                "score": score
            }

    # Clinical justification
    print("Clinical Rationale: The patient is in traumatic cardiac arrest, which is the most immediate life-threat.")
    print("This requires the highest level of surgical capability (Level 1 Trauma Center) for any chance of survival, even with a longer transport time.")
    print("The toxicology issue is secondary to the immediate need for surgical intervention.\n")

    # Output the final calculation as requested
    print(f"The best choice is the {best_hospital['name']}.")
    print("The final calculation for this choice is:")
    
    # Printing each number in the final equation
    trauma_component = best_hospital['trauma_cap'] * trauma_urgency
    tox_component = best_hospital['tox_cap'] * toxicology_urgency
    time_component = best_hospital['time'] * time_penalty_factor
    
    print(f"Score = (Trauma Capability[{best_hospital['trauma_cap']}] * Trauma Urgency[{trauma_urgency}]) + " \
          f"(Toxin Capability[{best_hospital['tox_cap']}] * Toxin Urgency[{toxicology_urgency}]) - " \
          f"(Transport Time[{best_hospital['time']}] * Time Penalty[{time_penalty_factor}])")

    print(f"Score = ({trauma_component}) + ({tox_component}) - ({time_component})")
    print(f"Final Score = {best_hospital['score']:.0f}")

solve_ems_scenario()
<<<E>>>