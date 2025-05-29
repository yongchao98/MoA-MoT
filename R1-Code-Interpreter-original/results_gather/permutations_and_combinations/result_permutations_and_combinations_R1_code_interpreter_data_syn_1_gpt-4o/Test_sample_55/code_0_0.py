# Define the positions of the books
positions = ["", "D", "C", "A", "G", "F", "", ""]

# Place B and E in the remaining positions
for i in range(1, 8):
    if positions[i] == "":
        positions[i] = "B"
        break

for i in range(1, 8):
    if positions[i] == "":
        positions[i] = "E"
        break

print(positions[1:])