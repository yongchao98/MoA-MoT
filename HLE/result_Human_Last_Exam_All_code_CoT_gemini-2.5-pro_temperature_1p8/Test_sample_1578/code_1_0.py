def oni_mirror_test():
    """
    This program explains the logical solution to the Mirror and the Oni puzzle.
    It demonstrates why using a second mirror is the most effective method.
    """
    print("Executing the logical test to expose the Oni (demon)...")
    print("-" * 50)
    
    print("The chosen strategy is: Use a second, trusted mirror to observe the suspected mirror.")
    print("The test procedure is as follows:")
    print("1. Turn your back to the suspected mirror (let's call it Mirror A).")
    print("2. In your hands, hold a second, smaller mirror (Mirror B).")
    print("3. Look into Mirror B to observe what is happening in Mirror A behind you.")
    print("-" * 50)
    
    print("Now, let's analyze the two possible outcomes:\n")
    
    # Case 1: The suspect is a genuine mirror
    print(">>> OUTCOME 1: If it is a REAL MIRROR <<<")
    print("A real mirror follows the laws of physics.")
    print("Since your back is facing Mirror A, Mirror A will reflect an image of your back.")
    print("When you look into Mirror B, you will see the reflection of Mirror A, which is showing your back.")
    print("RESULT: You see your own back. This is physically consistent and logical.\n")
    
    # Case 2: The suspect is a demon
    print(">>> OUTCOME 2: If it is a DEMON <<<")
    print("The demon must mimic your actions to maintain the illusion.")
    print("However, by turning your back, you have created a logical paradox for the demon:")
    print("  - CONTRADICTION A: To act like a proper reflection, the demon must show your back. But if it does this, it is NOT mimicking your CURRENT action (which is looking into Mirror B).")
    print("  - CONTRADICTION B: To mimic your CURRENT action, the demon would have to also pretend to look into a handheld mirror. But it cannot do this while simultaneously being a reflection of your back.")
    print("RESULT: The demon is caught. It cannot resolve the paradox. Its actions will not match a true reflection.")
    print("-" * 50)

# Run the logical demonstration
oni_mirror_test()