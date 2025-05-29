# Initialize the houses with empty dictionaries
houses = [{}, {}, {}, {}]

# Fill in the known information
houses[0]['name'] = 'Bob'
houses[0]['color'] = 'brown'
houses[0]['phone'] = 'Samsung Galaxy S21'

houses[1]['color'] = 'blue'
houses[1]['phone'] = 'Huawei P50'
houses[1]['name'] = 'Carol'

houses[2]['name'] = 'Alice'
houses[2]['phone'] = 'OnePlus 9'

# Deduce the remaining information
# Since Bob is in house 1, Carol is in house 2, and Alice is in house 3, Arnold must be in house 4
houses[3]['name'] = 'Arnold'

# From clue 8, the person who smokes Pall Mall is directly left of Alice, so they are in house 2
houses[1]['cigar'] = 'Pall Mall'

# From clue 3, the Prince smoker is directly left of the person who smokes Blue Master
# Since house 2 is Pall Mall, house 3 must be Blue Master, and house 2 must be Prince
houses[2]['cigar'] = 'Blue Master'
houses[1]['cigar'] = 'Prince'

# From clue 5, the Dunhill smoker loves white, so they must be in house 4
houses[3]['cigar'] = 'Dunhill'
houses[3]['color'] = 'white'

# The only remaining cigar is Pall Mall for house 1
houses[0]['cigar'] = 'Pall Mall'

# The only remaining phone is Google Pixel 6 for house 4
houses[3]['phone'] = 'Google Pixel 6'

# The only remaining color is purple for house 3
houses[2]['color'] = 'purple'

# Print the name of the person in house 1
print(houses[0]['name'])