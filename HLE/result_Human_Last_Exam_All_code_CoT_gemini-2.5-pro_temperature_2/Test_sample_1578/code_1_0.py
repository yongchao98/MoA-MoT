import random
import time

def detect_demon_in_mirror():
    """
    This function simulates the process of detecting a demon in a mirror
    by using rapid, unpredictable movements to expose its reaction lag.
    """
    
    actions = [
        "Raise Left Hand",
        "Duck Quickly",
        "Wink Right Eye",
        "Spin Around",
        "Look Up",
        "Stick Out Tongue"
    ]
    
    # Your current action and the action the demon is mimicking.
    # We start from a neutral state.
    your_action = "Stand Perfectly Still"
    demons_action = "Stand Perfectly Still"
    
    print("You are facing the mirror, trying to determine if it holds a demon.")
    print("You decide to test its reaction time with unpredictable movements.")
    print("-" * 60)
    
    # The test begins. We'll run for a few moments.
    for i in range(1, 6):
        print(f"Moment {i}:")
        
        # In a real scenario, the reflection would be instantaneous.
        # A demon, however, would still be copying your PREVIOUS move.
        # This simulates the demon's reaction lag.
        if i > 1: # On the first move, the demon will still be still
            demons_action = your_previous_action
            
        # You perform a new, random action.
        your_previous_action = your_action
        your_action = random.choice(actions)

        print(f"    - Your Action:          You '{your_action}'")
        print(f"    - 'Reflection's' Action: The reflection '{demons_action}'")
        
        # Compare your current action to the reflection's action.
        if your_action != demons_action:
            print("\n>>> A-HA! A MISMATCH! <<<")
            print(f"Your move was '{your_action}', but the reflection performed '{demons_action}'.")
            print("The lag in its reaction proves this is no reflection, but a demon's illusion!")
            print("-" * 60)
            return
            
        print("    The actions match... for now. You continue the test.")
        print("-" * 60)
        time.sleep(1) # Pause to make the output readable.

detect_demon_in_mirror()