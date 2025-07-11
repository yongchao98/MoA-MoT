import time
import random

def test_the_mirror():
    """
    This program simulates a test to determine if a reflection is real
    or a mimicry by a demon (Oni) with a reaction delay.
    It demonstrates the principle of exploiting latency in a mimic.
    """

    movements = [
        "Raise right hand",
        "Wink left eye",
        "Snap fingers",
        "Turn head sharply right",
        "Look down quickly",
        "Puff out cheeks"
    ]
    # Let's assume the previous state was being still.
    previous_action = "Standing still"

    print("You suspect the mirror houses a demon. You decide to test it with fast, erratic movements.")
    print("A true reflection is bound by the laws of physics (instantaneous).")
    print("A demon is bound by perception and reaction, which takes time.")
    print("-" * 70)

    # The simulation runs for 3 quick movements.
    for test_number in range(1, 4):
        # You perform a new action at Time T=0.
        current_action = random.choice(movements)
        print(f"\nMOVEMENT #{test_number} (Time: T=0.0s)")
        print(f"Your real action: '{current_action}'")

        # --- Check 1: The True Mirror ---
        # The reflection in a true mirror is instantaneous.
        true_reflection = current_action
        print("\nIf it is a True Mirror, the 'equation' of your action and its reflection is:")
        # In this logical equation, we output the test number.
        print(f"  Test {test_number}: Does Your Action ('{current_action}') == True Reflection ('{true_reflection}')?")
        print(f"  Result: {current_action == true_reflection}. It is a perfect, instant match.")


        # --- Check 2: The Demon Mirror ---
        # The demon is still processing your *last* action at the moment you begin your new one.
        demon_reflection = previous_action
        print("\nIf it is a Demon Mirror, the 'equation' is:")
        # In this logical equation, we also output the test number.
        print(f"  Test {test_number}: Does Your Action ('{current_action}') == Demon Reflection ('{demon_reflection}')?")
        print(f"  Result: {current_action == demon_reflection}. A clear mismatch is detected!")

        # Update the state for the next iteration of the loop.
        previous_action = current_action
        print("-" * 70)

    print("\nCONCLUSION: The consistent lag between your actions and what you see proves")
    print("that the figure in the mirror is not a reflection, but a demonic mimic.")

# Run the simulation.
test_the_mirror()