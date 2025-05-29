# Define the houses and their attributes
houses = [None] * 4

# Assign known values based on the clues
houses[1] = {'name': 'Arnold', 'nationality': 'brit', 'lunch': 'soup'}
houses[2] = {'name': 'Bob', 'nationality': 'german', 'child': 'Bella', 'lunch': 'stir fry'}
houses[0] = {'name': 'Carol', 'nationality': 'chinese', 'child': 'Alice'}

# Determine the remaining attributes
for i in range(4):
    if houses[i] is None:
        houses[i] = {'name': 'alice', 'nationality': 'dane', 'child': 'Billy', 'lunch': 'pizza'}

# Find the name of the person in House 1
house_1_name = houses[0]['name']
print(house_1_name)