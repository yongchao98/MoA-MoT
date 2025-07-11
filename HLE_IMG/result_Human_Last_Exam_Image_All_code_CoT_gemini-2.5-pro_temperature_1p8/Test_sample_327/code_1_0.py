# Step 1: Define the properties for each nanoribbon based on visual analysis.

# Ribbon 1 properties: Armchair, Width=8, Semiconducting
edge1 = 'A'
width1 = 8
band1 = 1

# Ribbon 2 properties: Armchair, Width=7, Metallic
edge2 = 'A'
width2 = 7
band2 = 0

# Ribbon 3 properties: Armchair, Width=6, Semiconducting
edge3 = 'A'
width3 = 6
band3 = 1

# Step 2: Create the classification string for each ribbon.
# The format is "Edge_Width_Band".
class1 = f"{edge1}{width1}{band1}"
class2 = f"{edge2}{width2}{band2}"
class3 = f"{edge3}{width3}{band3}"

# Step 3: Concatenate the individual classifications to get the final result.
# The code shows how each part, including each number, contributes to the final string.
final_result = class1 + class2 + class3

# Step 4: Print the final concatenated string.
print(final_result)