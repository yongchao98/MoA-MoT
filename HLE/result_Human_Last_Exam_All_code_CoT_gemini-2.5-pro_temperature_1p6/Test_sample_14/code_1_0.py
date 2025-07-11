def find_stunning_t1_flyer():
    """
    Identifies and prints the T1 flying unit in BAR that can stun enemies.
    """
    # Define the characteristics of the unit based on game knowledge.
    unit_tier = 1
    unit_name = "Hawk"
    unit_faction = "ARM"
    unit_ability = "Stun/Paralyze"

    print(f"The flying unit from a Tier {unit_tier} building in BAR that can shoot and stun enemies is the {unit_faction} {unit_name}.")
    print("\nThis unit is a gunship known for its ability to paralyze enemy units with its laser weapon, making them unable to move or fire.")

    # Fulfilling the request to show each part of the "final equation".
    print("\n--- Answer Breakdown ---")
    print(f'Tier({unit_tier}) + Unit("{unit_name}") = Ability("{unit_ability}")')

# Execute the function to provide the answer.
find_stunning_t1_flyer()