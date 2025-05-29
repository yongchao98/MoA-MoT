# Initialize the houses with empty dictionaries
houses = [{} for _ in range(4)]

# Apply the clues
# Clue 1: The person's child named Alice is in the second house.
houses[1]['child'] = 'alice'

# Clue 8: Arnold's child is named Bella.
# This means Arnold cannot be in house 2, as house 2's child is Alice.
# Arnold's child is Bella, so Arnold is not in house 2.

# Clue 6: Carol is directly left of Arnold.
# This means Arnold cannot be in house 1, as there is no house to the left of house 1.
# Therefore, Arnold must be in house 3 or 4, and Carol must be in house 2 or 3.

# Clue 5: Bob is directly left of the person who loves the bouquet of iris.
# This means Bob cannot be in house 4, as there is no house to the right of house 4.
# Therefore, Bob must be in house 1, 2, or 3.

# Clue 7: The person who loves the bouquet of iris is the person who smokes Blue Master.
# This means the person who loves iris cannot be in house 1, as Bob is in house 1 and is directly left of the person who loves iris.

# Clue 9: Alice and the person who loves the bouquet of lilies are next to each other.
# This means the person who loves lilies cannot be in house 1, as Alice is in house 2.

# Let's try to deduce the positions based on these constraints.
# We know:
# - House 2's child is Alice.
# - Arnold is in house 3 or 4.
# - Carol is in house 2 or 3.
# - Bob is in house 1, 2, or 3.
# - The person who loves iris is in house 2, 3, or 4.
# - The person who loves lilies is in house 2 or 3.

# Let's try to place Arnold and Carol first.
# If Arnold is in house 3, then Carol must be in house 2.
# If Arnold is in house 4, then Carol must be in house 3.

# Let's try Arnold in house 3 and Carol in house 2.
houses[2]['name'] = 'carol'
houses[3]['name'] = 'arnold'
houses[3]['child'] = 'bella'

# Now, Bob must be in house 1 or 2.
# If Bob is in house 1, then the person who loves iris is in house 2.
# But house 2 is Carol, and we don't have enough information to place iris there.
# Let's try Bob in house 1.
houses[0]['name'] = 'bob'

# Now, the person who loves iris must be in house 2 or 3.
# If the person who loves iris is in house 3, then Bob is directly left of them, which is correct.
# Let's place iris in house 3.
houses[3]['flower'] = 'iris'
houses[3]['cigar'] = 'blue master'

# Now, we have:
# House 1: Bob
# House 2: Carol, child Alice
# House 3: Arnold, child Bella, flower iris, cigar blue master

# The person who loves lilies must be in house 2, as they are next to Alice.
houses[2]['flower'] = 'lilies'

# The remaining person is Alice, who must be in house 4.
houses[1]['name'] = 'alice'

# Now, we have:
# House 1: Bob
# House 2: Carol, child Alice, flower lilies
# House 3: Arnold, child Bella, flower iris, cigar blue master
# House 4: Alice

# The remaining cigars are Prince, Dunhill, and Pall Mall.
# The remaining children are Billy and Timothy.
# The remaining flowers are carnations and daffodils.

# Clue 2: The person who is the mother of Billy and the Prince smoker are next to each other.
# Clue 3: The person who is the mother of Timothy and the Dunhill smoker are next to each other.
# Clue 4: The person who loves a bouquet of daffodils is the mother of Billy.

# Since Arnold is in house 3 and smokes Blue Master, the Prince smoker must be in house 2 or 4.
# Since Carol is in house 2, the Prince smoker must be in house 4.
houses[3]['cigar'] = 'prince'

# The person who loves daffodils is the mother of Billy, and they must be next to the Prince smoker.
# Since the Prince smoker is in house 4, the person who loves daffodils must be in house 3.
houses[2]['flower'] = 'daffodils'
houses[2]['child'] = 'billy'

# The remaining cigar is Dunhill, which must be in house 2.
houses[1]['cigar'] = 'dunhill'

# The remaining child is Timothy, who must be in house 1.
houses[0]['child'] = 'timothy'

# The remaining flower is carnations, which must be in house 1.
houses[0]['flower'] = 'carnations'

# Now, we have:
# House 1: Bob, child Timothy, flower carnations
# House 2: Carol, child Billy, flower daffodils, cigar Dunhill
# House 3: Arnold, child Bella, flower iris, cigar Blue Master
# House 4: Alice, cigar Prince

# The name of the person who lives in House 1 is Bob.
print(houses[0]['name'])