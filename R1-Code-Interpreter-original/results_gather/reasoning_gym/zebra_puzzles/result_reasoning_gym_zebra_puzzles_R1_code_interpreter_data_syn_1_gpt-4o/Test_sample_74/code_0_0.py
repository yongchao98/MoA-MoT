# Initialize the houses with empty dictionaries
houses = [{'name': None, 'cigar': None, 'lunch': None, 'drink': None} for _ in range(4)]

# Apply the clues
# Clue 7: The person partial to Pall Mall is in the second house.
houses[1]['cigar'] = 'pall mall'

# Clue 8: The coffee drinker is Bob.
for house in houses:
    if house['name'] is None:
        house['name'] = 'bob'
        house['drink'] = 'coffee'
        break

# Clue 3: Bob is directly left of Arnold.
for i in range(3):
    if houses[i]['name'] == 'bob':
        houses[i+1]['name'] = 'arnold'
        break

# Clue 6: The Prince smoker is the person who likes milk.
for house in houses:
    if house['cigar'] is None and house['drink'] is None:
        house['cigar'] = 'prince'
        house['drink'] = 'milk'
        break

# Clue 1: The person who is a pizza lover is the one who only drinks water.
for house in houses:
    if house['lunch'] is None and house['drink'] is None:
        house['lunch'] = 'pizza'
        house['drink'] = 'water'
        break

# Clue 2: The person who smokes Blue Master is the person who loves stir fry.
for house in houses:
    if house['cigar'] is None and house['lunch'] is None:
        house['cigar'] = 'blue master'
        house['lunch'] = 'stir fry'
        break

# Clue 9: The person who loves stir fry and Carol are next to each other.
for i in range(3):
    if houses[i]['lunch'] == 'stir fry' and houses[i+1]['name'] is None:
        houses[i+1]['name'] = 'carol'
        break
    elif houses[i+1]['lunch'] == 'stir fry' and houses[i]['name'] is None:
        houses[i]['name'] = 'carol'
        break

# Clue 4: The Dunhill smoker and the person who loves the soup are next to each other.
# Clue 5: The person who loves the soup and the boba tea drinker are next to each other.
for i in range(3):
    if houses[i]['cigar'] is None and houses[i+1]['lunch'] is None:
        houses[i]['cigar'] = 'dunhill'
        houses[i+1]['lunch'] = 'soup'
        houses[i+1]['drink'] = 'boba tea'
        break
    elif houses[i+1]['cigar'] is None and houses[i]['lunch'] is None:
        houses[i+1]['cigar'] = 'dunhill'
        houses[i]['lunch'] = 'soup'
        houses[i]['drink'] = 'boba tea'
        break

# Fill in the remaining name
for house in houses:
    if house['name'] is None:
        house['name'] = 'alice'

# Determine the name of the person in House 1
house_1_name = houses[0]['name']
print(house_1_name)