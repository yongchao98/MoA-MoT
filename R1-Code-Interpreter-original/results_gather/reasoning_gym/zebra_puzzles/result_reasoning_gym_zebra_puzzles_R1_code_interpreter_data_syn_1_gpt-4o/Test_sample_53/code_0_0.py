# Initialize the houses with empty dictionaries
houses = [{}, {}, {}, {}]

# Place known information
houses[2]['cigar'] = 'Pall Mall'  # Clue 7
houses[2]['color'] = 'blue'  # Clue 9
houses[2]['name'] = 'Arnold'  # Clue 9

# Alice uses Google Pixel 6 and smokes Blue Master
for house in houses:
    if 'name' not in house:
        house['name'] = 'Alice'
        house['phone'] = 'Google Pixel 6'
        house['cigar'] = 'Blue Master'
        break

# The person who loves brown uses OnePlus 9
for i, house in enumerate(houses):
    if 'color' not in house and 'phone' not in house:
        house['color'] = 'brown'
        house['phone'] = 'OnePlus 9'
        if i > 0:
            houses[i-1]['color'] = 'white'  # Clue 8
        break

# The Prince smoker uses Samsung Galaxy S21 and loves purple
for house in houses:
    if 'cigar' not in house:
        house['cigar'] = 'Prince'
        house['phone'] = 'Samsung Galaxy S21'
        house['color'] = 'purple'
        break

# Carol is next to the person who loves brown
for i, house in enumerate(houses):
    if house.get('color') == 'brown':
        if i > 0 and 'name' not in houses[i-1]:
            houses[i-1]['name'] = 'Carol'
        elif i < 3 and 'name' not in houses[i+1]:
            houses[i+1]['name'] = 'Carol'
        break

# The remaining person is Bob
for house in houses:
    if 'name' not in house:
        house['name'] = 'Bob'

# Find the name of the person in House 1
print(houses[0]['name'])