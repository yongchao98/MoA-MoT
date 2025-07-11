def solve_magic_puzzle():
    """
    Calculates the maximum life Player B will lose based on the described
    Magic: The Gathering board state and optimal play from both sides.
    """

    print("To maximize damage, Player A will perform a series of actions before and during combat.")
    print("--- Pre-Combat Main Phase ---")
    print("1. Player A casts 'Iron Apprentice' (1/1) for {1} mana.")
    print("2. 'Replication Specialist' triggers, and Player A pays {2} mana to create a copy, another 1/1 'Iron Apprentice'.")
    
    print("\n--- Combat Phase ---")
    print("3. Player A crews 'Mukotai Soulripper' with a 2/2 Ninja token. The Ninja cannot attack.")
    print("4. Player A attacks with all other creatures.")
    print("5. On attack, Player A sacrifices an 'Iron Apprentice' to 'Mukotai Soulripper', making it a 4/5 with menace.")
    print("6. The sacrificed apprentice's ability moves its +1/+1 counter to the other apprentice, making it a 2/2.")

    # The final set of creatures attacking and their powers.
    attackers = {
        "Mukotai Soulripper (menace)": 4,
        "Ironhoof Boar": 4,
        "Replication Specialist": 2,
        "Scrap Welder (enchanted)": 2,
        "Ninja Token": 2,
        "Iron Apprentice (upgraded)": 2,
    }
    total_attacking_power = sum(attackers.values())
    print(f"\nPlayer A is attacking with a total of {total_attacking_power} power.")

    print("\n--- Declare Blockers Phase ---")
    print("Player B has two 2/1 blockers and must minimize life loss.")
    print("Player B's optimal play is to block the two highest-power non-menace creatures.")
    
    # Player B's blocking decisions
    damage_prevented_boar = 4
    damage_prevented_other = 2
    print(f"7. Player B blocks 'Ironhoof Boar', preventing {damage_prevented_boar} damage.")
    print(f"8. Player B blocks one of the 2-power creatures (e.g., Scrap Welder), preventing {damage_prevented_other} damage.")
    
    # Unblocked attackers deal damage.
    unblocked_attackers_damage = {
        "Mukotai Soulripper": 4,
        "Replication Specialist": 2,
        "Ninja Token": 2,
        "Iron Apprentice": 2,
    }
    
    # Calculate and display the final damage
    total_life_loss = sum(unblocked_attackers_damage.values())
    
    print("\n--- Damage Calculation ---")
    print("The unblocked creatures deal their damage to Player B.")

    # Build the final equation string dynamically
    damage_values = list(unblocked_attackers_damage.values())
    equation_parts = [str(val) for val in damage_values]
    equation_string = " + ".join(equation_parts)
    
    print(f"Total life loss = {equation_string} = {total_life_loss}")

solve_magic_puzzle()
<<<10>>>