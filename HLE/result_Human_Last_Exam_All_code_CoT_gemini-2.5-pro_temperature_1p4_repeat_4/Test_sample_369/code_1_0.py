# The line of Latin hexameter is:
# "verum ubi equi atque hominis casu convenit imago."
# The scansion represents each of the six feet as 'D' for a Dactyl or 'S' for a Spondee.

# Define the six feet of the scanned line
foot1 = "D"
foot2 = "S"
foot3 = "S"
foot4 = "S"
foot5 = "D"
foot6 = "S"

# Create a list of the feet
scansion = [foot1, foot2, foot3, foot4, foot5, foot6]

# Print the final scansion, showing each foot in the sequence
print("The scansion of the line is:")
print(f"Foot 1: {scansion[0]}")
print(f"Foot 2: {scansion[1]}")
print(f"Foot 3: {scansion[2]}")
print(f"Foot 4: {scansion[3]}")
print(f"Foot 5: {scansion[4]}")
print(f"Foot 6: {scansion[5]}")

print("\nFinal Result:")
# Use the join function to create the final equation format
final_result = " + ".join(scansion)
print(final_result)