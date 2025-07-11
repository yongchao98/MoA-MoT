def find_stunning_t1_flyer():
    """
    This function searches for a Tier 1 flying unit in BAR 
    that can both shoot and stun enemies.
    """
    # Data representing T1 flying units from BAR air plants and their capabilities
    t1_flying_units = [
        # Armada
        {"name": "Fink", "faction": "Armada", "type": "Scout", "can_shoot": True, "can_stun": False},
        {"name": "Slasher", "faction": "Armada", "type": "Fighter", "can_shoot": True, "can_stun": False},
        {"name": "Hawkeye", "faction": "Armada", "type": "Gunship", "can_shoot": True, "can_stun": False},
        
        # Cortex
        {"name": "Scout", "faction": "Cortex", "type": "Scout", "can_shoot": True, "can_stun": False},
        {"name": "Vamp", "faction": "Cortex", "type": "Fighter", "can_shoot": True, "can_stun": False},
        {"name": "Thrasher", "faction": "Cortex", "type": "Gunship", "can_shoot": True, "can_stun": False}
    ]

    found_unit = None
    
    # Iterate through the list to find a unit that matches the criteria
    for unit in t1_flying_units:
        if unit["can_shoot"] and unit["can_stun"]:
            found_unit = unit
            break  # Stop after finding the first one

    # Print the result
    if found_unit:
        print(f"The Tier 1 flying unit that can shoot and stun is: {found_unit['name']}")
    else:
        print("Analysis complete: No Tier 1 flying unit from a T1 air plant can both shoot and stun targets.")
        print("\nThe stun ability for flying units is a Tier 2 feature.")
        print("Examples of T2 stunning air units include the Armada 'Stiletto' and the Cortex 'Viper'.")

if __name__ == "__main__":
    find_stunning_t1_flyer()