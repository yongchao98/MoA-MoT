import random

def solve_mirror_puzzle():
    """
    Simulates a test to distinguish a true reflection from an Oni by measuring reaction time lag.
    """
    # Step 1: Define the timing of the events.
    # The user performs a sudden, unpredictable action at time t=0.
    user_action_time = 0.0

    # A true reflection appears almost instantly. The time is based on the speed of light
    # over a short distance, which is negligible for human perception.
    light_travel_time = 0.0

    # The Oni, however fast, has a biological or magical reaction time.
    # Let's assume a superhuman, but still non-zero, reaction time of 50 milliseconds.
    oni_reaction_time = 0.05  # in seconds

    # Step 2: Calculate when the user would see each event.
    reflection_action_time = user_action_time + light_travel_time
    oni_action_time = user_action_time + oni_reaction_time

    # Step 3: Use the final equation to calculate the detectable lag.
    lag = oni_action_time - reflection_action_time

    # Step 4: Print the analysis and the results.
    print("Testing for the Oni in the mirror by checking for reaction time.")
    print("A quick, unpredictable movement is performed at Time = 0.000s.\n")

    print("--- SCENARIO 1: A True Reflection ---")
    print(f"A real reflection is instantaneous (travels at the speed of light).")
    print(f"Reflection appears at: {reflection_action_time:.3f}s")
    print("-" * 40)

    print("--- SCENARIO 2: An Oni's Imitation ---")
    print(f"The Oni must perceive the action and then react, causing a delay.")
    print(f"Assumed Oni Reaction Time: {oni_reaction_time:.3f}s")
    print(f"Oni imitation appears at: {oni_action_time:.3f}s")
    print("-" * 40)

    print("--- CONCLUSION: THE LAG EQUATION ---")
    print("The final equation to detect the Oni is: Lag = Oni's Action Time - Reflection's Action Time")
    print("The numbers in the final equation are:")
    # Printing each number in the final equation as requested.
    print(f"{lag:.3f}s = {oni_action_time:.3f}s - {reflection_action_time:.3f}s")

    if lag > 0:
        print("\nResult: A lag was detected! The figure is an illusion created by the Oni.")
    else:
        print("\nResult: No lag was detected. It is a true reflection.")

# Run the simulation
solve_mirror_puzzle()