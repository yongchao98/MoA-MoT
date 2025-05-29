# Initialize the houses with empty dictionaries
houses = [{}, {}, {}, {}]

# Fill in the known information
houses[0]['name'] = 'carol'
houses[0]['cigar'] = 'prince'
houses[0]['lunch'] = 'soup'
houses[0]['color'] = 'brown'

houses[2]['lunch'] = 'pizza'

houses[1]['name'] = 'bob'
houses[1]['lunch'] = 'stir fry'
houses[1]['color'] = 'blue'

houses[3]['name'] = 'arnold'
houses[3]['cigar'] = 'dunhill'

# Deduce the remaining information
# House 1 is Carol, so House 2 must be the one with grilled cheese and purple
houses[1]['lunch'] = 'grilled cheese'
houses[1]['color'] = 'purple'

# House 2 is directly left of the Pall Mall smoker, so House 3 must be Pall Mall
houses[2]['cigar'] = 'pall mall'

# The only remaining cigar for House 4 is Blue Master
houses[3]['cigar'] = 'blue master'

# The only remaining name for House 2 is Alice
houses[1]['name'] = 'alice'

# The only remaining lunch for House 4 is pizza
houses[3]['lunch'] = 'pizza'

# The only remaining color for House 3 is white
houses[2]['color'] = 'white'

# Print the name of the person in House 1
print(houses[0]['name'])