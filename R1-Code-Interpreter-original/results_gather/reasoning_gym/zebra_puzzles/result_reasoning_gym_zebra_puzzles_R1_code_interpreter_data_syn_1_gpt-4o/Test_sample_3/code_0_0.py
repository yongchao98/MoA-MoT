# Let's create a list of houses with their attributes
houses = [None] * 4

# We know:
# - Bob loves white and smokes Prince
# - Alice loves brown and smokes Blue Master
# - The coffee drinker is in the third house
# - The person who loves blue drinks water and smokes Pall Mall

# Let's fill in the information we know:
houses[2] = {'drink': 'coffee'}  # Third house drinks coffee

# Bob loves white and smokes Prince
# Alice loves brown and smokes Blue Master
# The person who loves blue drinks water and smokes Pall Mall

# Let's try to place Bob and Alice
for i in range(4):
    if houses[i] is None:
        houses[i] = {'name': 'Bob', 'color': 'white', 'cigar': 'Prince'}
        break

for i in range(4):
    if houses[i] is None:
        houses[i] = {'name': 'Alice', 'color': 'brown', 'cigar': 'Blue Master'}
        break

# Now, let's place the person who loves blue
for i in range(4):
    if houses[i] is None:
        houses[i] = {'color': 'blue', 'drink': 'water', 'cigar': 'Pall Mall'}
        break

# The remaining house must be the one with purple and milk
for i in range(4):
    if houses[i] is None:
        houses[i] = {'color': 'purple', 'drink': 'milk'}

# Now, let's determine the names based on the remaining clues
# The person who loves purple and Arnold are next to each other
# The person who loves purple drinks milk

# Let's find the house with purple
purple_index = next(i for i, house in enumerate(houses) if house['color'] == 'purple')

# Arnold must be next to the person who loves purple
if purple_index > 0 and 'name' not in houses[purple_index - 1]:
    houses[purple_index - 1]['name'] = 'Arnold'
elif purple_index < 3 and 'name' not in houses[purple_index + 1]:
    houses[purple_index + 1]['name'] = 'Arnold'

# The remaining name is Carol
for house in houses:
    if 'name' not in house:
        house['name'] = 'Carol'

# Find the name of the person in House 1
house_1_name = houses[0]['name']
print(house_1_name)