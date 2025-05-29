# Initialize houses with all possible keys
houses = [{'name': None, 'color': None, 'drink': None, 'cigar': None} for _ in range(4)]

# Fill in the known values
houses[2]['drink'] = 'coffee'  # Third house drinks coffee

# Bob loves white and smokes Prince
# Alice loves brown and smokes Blue Master
# The person who loves blue drinks water and smokes Pall Mall

# Place Bob and Alice
for i in range(4):
    if houses[i]['name'] is None:
        houses[i].update({'name': 'Bob', 'color': 'white', 'cigar': 'Prince'})
        break

for i in range(4):
    if houses[i]['name'] is None:
        houses[i].update({'name': 'Alice', 'color': 'brown', 'cigar': 'Blue Master'})
        break

# Place the person who loves blue
for i in range(4):
    if houses[i]['color'] is None:
        houses[i].update({'color': 'blue', 'drink': 'water', 'cigar': 'Pall Mall'})
        break

# The remaining house must be the one with purple and milk
for i in range(4):
    if houses[i]['color'] is None:
        houses[i].update({'color': 'purple', 'drink': 'milk'})

# Determine the names based on the remaining clues
# The person who loves purple and Arnold are next to each other
# The person who loves purple drinks milk

# Find the house with purple
purple_index = next(i for i, house in enumerate(houses) if house['color'] == 'purple')

# Arnold must be next to the person who loves purple
if purple_index > 0 and houses[purple_index - 1]['name'] is None:
    houses[purple_index - 1]['name'] = 'Arnold'
elif purple_index < 3 and houses[purple_index + 1]['name'] is None:
    houses[purple_index + 1]['name'] = 'Arnold'

# The remaining name is Carol
for house in houses:
    if house['name'] is None:
        house['name'] = 'Carol'

# Find the name of the person in House 1
house_1_name = houses[0]['name']
print(house_1_name)