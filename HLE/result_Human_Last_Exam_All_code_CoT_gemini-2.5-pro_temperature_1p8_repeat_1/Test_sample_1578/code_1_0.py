def solve_mirror_puzzle():
    """
    This program outlines a logical procedure to solve the Japanese puzzle
    of the Mirror and the Oni, determining if a reflection is real or a demon's illusion.
    """
    print("The Puzzle: You are looking into a mirror. How do you determine if the figure you see is your reflection or an Oni (demon) perfectly mimicking you?")
    print("-" * 80)
    print("The Chosen Method (G): Use an opaque barrier to test the nature of the reflection.")
    print("This test is effective because a true reflection is bound by the laws of optics, while a demon mimic is bound by its own perception.")
    print("-" * 80)
    
    print("\nStep 1: Stand in front of the alleged mirror so you can see your 'reflection' clearly.")
    
    print("\nStep 2: Procure a large, opaque barrier (e.g., a sturdy piece of cardboard or a thick book).")

    print("\nStep 3: Hold the barrier up to partially obscure your body from the mirror. For instance, block the view of your torso and your right hand, but ensure you can still see the figure's head and the figure can see yours.")
    
    print("\nStep 4: While your right hand is completely hidden behind the barrier (from the mirror's perspective), perform a specific, simple action. For example, raise three fingers.")

    print("\nStep 5: Carefully observe the 'reflection' of your hand. Since it's supposed to be obscured by the barrier in the reflection as well, you must deduce its action by what is not happening. The demon may give itself away.")
    
    print("\n--- ANALYSIS OF OUTCOMES ---")
    
    print("\nPossibility #1: The figure in the mirror perfectly mimics the obscuring of your hand. You see the barrier in the mirror, and the hand behind it is not visible. Nothing happens to give away a hidden movement.")
    print(">>> Conclusion: This is behaving exactly as a true reflection would. The mirror cannot show what it cannot 'see'. The object is most likely a real mirror.")
    
    print("\nPossibility #2: The figure in the mirror makes a mistake. For example, its hand remains still and visible when it should be obscured, or it fails to perform your hidden action, revealing a static hand if you peek around the barrier.")
    print(">>> Conclusion: You have found the demon! It cannot mimic an action it cannot see. Its illusion is flawed because it relies on observing you, and you have successfully blocked its line of sight. It is not a reflection.")

# Execute the logical walkthrough
solve_mirror_puzzle()
