def find_stunning_t1_flyer():
    """
    Finds and prints the T1 flying unit in BAR that can shoot and stun.
    """
    # Dataset of T1 flying units from BAR's air plants
    t1_air_units = [
        # Armada
        {"name": "Flea", "faction": "Armada", "can_shoot": True, "can_stun": False},
        {"name": "Fink", "faction": "Armada", "can_shoot": True, "can_stun": False},
        {"name": "Hawkeye", "faction": "Armada", "can_shoot": True, "can_stun": False},
        
        # Cortex
        {"name": "Vulture", "faction": "Cortex", "can_shoot": True, "can_stun": False},
        {"name": "Thrasher", "faction": "Cortex", "can_shoot": True, "can_stun": True},
        {"name": "Shuriken", "faction": "Cortex", "can_shoot": True, "can_stun": False}
    ]

    print("Searching for a Tier 1 flying unit that can shoot and stun...")
    
    found_unit = False
    for unit in t1_air_units:
        # Check if the unit meets the criteria
        if unit["can_shoot"] and unit["can_stun"]:
            print(f"Found unit: {unit['name']}")
            found_unit = True
            
    if not found_unit:
        print("No Tier 1 flying unit found that meets the criteria.")

find_stunning_t1_flyer()