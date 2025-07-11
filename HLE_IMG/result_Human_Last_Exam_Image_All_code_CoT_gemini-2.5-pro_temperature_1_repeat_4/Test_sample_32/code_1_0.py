def explain_moss_difference():
    """
    This function explains the key difference between the sporophytes of
    Grimmia montana and Grimmia curviseta.
    """
    
    # Define the characteristics for each species
    g_montana = {
        "name": "Grimmia montana",
        "seta": "very short and straight",
        "capsule_position": "immersed among the leaves"
    }

    g_curviseta = {
        "name": "Grimmia curviseta",
        "seta": "longer and strongly curved (cygneous or swan-necked)",
        "capsule_position": "emergent or exserted from the leaves"
    }

    # Print the explanation
    print("The main difference in the sporophyte that helps distinguish Grimmia montana and Grimmia curviseta at a glance is the seta (the stalk supporting the capsule):\n")
    
    print(f"1. In {g_montana['name']}:")
    print(f"   - The seta is {g_montana['seta']}.")
    print(f"   - This makes the capsule appear {g_montana['capsule_position']}.\n")
    
    print(f"2. In {g_curviseta['name']}:")
    print(f"   - The seta is {g_curviseta['seta']}.")
    print(f"   - This causes the capsule to be {g_curviseta['capsule_position']}.\n")

    print("Therefore, the most obvious distinguishing feature is the straight, short seta in G. montana versus the prominent, curved seta in G. curviseta.")

# Run the function to display the answer
explain_moss_difference()