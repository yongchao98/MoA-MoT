def find_bar_t1_stunning_flyers():
    """
    Identifies and prints the T1 flying units in the game BAR
    that have a weapon that can stun enemies.
    """
    # Data for the T1 flying units that can shoot and stun.
    # These are the gunships for each faction.
    units = [
        {
            "name": "Owl",
            "faction": "Armada",
            "building": "T1 Aeroplant",
            "description": "A gunship that fires a continuous laser beam, which stuns and damages light ground units."
        },
        {
            "name": "Locust",
            "faction": "Cortex",
            "building": "T1 Air Factory",
            "description": "A gunship that uses a stunning laser weapon, effective against ground targets."
        }
    ]

    print("The T1 flying units in BAR that can shoot and stun are:")
    for unit in units:
        print(f"\nUnit: {unit['name']}")
        print(f"Faction: {unit['faction']}")
        print(f"Produced from: {unit['building']}")
        print(f"Description: {unit['description']}")

# Execute the function to print the information.
find_bar_t1_stunning_flyers()