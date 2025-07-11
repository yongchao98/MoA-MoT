import time

def solve_mirror_puzzle():
    """
    This script simulates the test to determine if the figure in a mirror
    is a reflection or a demon, based on the principle of information asymmetry.
    It demonstrates the logic behind option G.
    """

    # --- Initial State ---
    # The 'equation' of our logical proof is built on these states:
    # State 1: Visible action (what the demon can see)
    # State 2: Hidden action (what the demon cannot see)
    # State 3: The expected action of a true reflection
    # State 4: The demon's mimicked action
    # The proof is complete when State 3 != State 4.
    
    person_visible_action = "Standing still, holding a large opaque barrier."
    person_hidden_action_initial = "None"
    
    print("--- The Oni in the Mirror Test ---")
    print("A simulation to determine if your 'reflection' is an imposter.")
    print("\n[STEP 1: The Setup]")
    print(f"You stand before the mirror. Your visible action (State 1) is: '{person_visible_action}'")
    print(f"Your hidden action behind the barrier (State 2) is: '{person_hidden_action_initial}'")

    # The True Reflection's hidden state will perfectly mirror your hidden state (with L/R reversal).
    true_reflection_hidden_action = "None"
    # The Demon can only see your visible state, so it cannot know your hidden action.
    demon_hidden_action = "None"

    print(f"The True Reflection's hidden action (State 3) is: '{true_reflection_hidden_action}'")
    print(f"The Demon's mimicked hidden action (State 4) is: '{demon_hidden_action}'")
    print("RESULT: Currently, State 3 equals State 4. The demon is safe.")
    print("-" * 25)
    time.sleep(2)

    # --- The Test ---
    print("\n[STEP 2: The Hidden Action]")
    person_hidden_action_new = "Waving your Right Hand"
    print(f"Behind the barrier, you change your hidden action (State 2) to: '{person_hidden_action_new}'")
    
    # Update the states based on the new hidden action.
    # A true reflection would mirror this perfectly.
    true_reflection_hidden_action = "Waving its Left Hand"
    # The demon CANNOT see this change. Its knowledge is limited to State 1, which hasn't changed.
    demon_hidden_action = "None" # The demon continues to do nothing as it sees no change.

    print(f"The True Reflection's hidden action (State 3) is now: '{true_reflection_hidden_action}'")
    print(f"The Demon's mimicked hidden action (State 4) is still: '{demon_hidden_action}'")
    print("RESULT: Now, State 3 does NOT equal State 4.")
    print("-" * 25)
    time.sleep(2)

    # --- The Conclusion ---
    print("\n[STEP 3: The Revelation]")
    # Here we evaluate the final 'equation' of our logic test.
    print(f"Final logical check: Is ('{true_reflection_hidden_action}') != ('{demon_hidden_action}')?")
    if true_reflection_hidden_action != demon_hidden_action:
        print("YES. The actions do not match!")
        print("The figure in the mirror failed to replicate your hidden movement. It's a demon!")
    else:
        # This case shouldn't be reached in this simulation.
        print("The actions match. The nature of the figure is still unknown.")

solve_mirror_puzzle()
print("\n<<<G>>>")