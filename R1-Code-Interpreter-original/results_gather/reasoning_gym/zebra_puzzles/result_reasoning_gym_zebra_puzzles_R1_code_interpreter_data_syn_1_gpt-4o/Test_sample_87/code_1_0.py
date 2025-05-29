from itertools import permutations

# Define the possible options for each category
names = ['alice', 'arnold', 'carol', 'bob']
children = ['timothy', 'bella', 'billy', 'alice']
smoothies = ['desert', 'dragonfruit', 'darkness', 'butterscotch']
lunches = ['stir fry', 'grilled cheese', 'soup', 'pizza']

# Iterate over all permutations of the names
for name_perm in permutations(names):
    # Assign names to houses
    house1, house2, house3, house4 = name_perm
    
    # Check if Bob is the mother of Timothy
    if house1 == 'bob' and children[0] != 'timothy':
        continue
    if house2 == 'bob' and children[1] != 'timothy':
        continue
    if house3 == 'bob' and children[2] != 'timothy':
        continue
    if house4 == 'bob' and children[3] != 'timothy':
        continue
    
    # Check if Arnold loves the soup
    if house1 == 'arnold' and lunches[0] != 'soup':
        continue
    if house2 == 'arnold' and lunches[1] != 'soup':
        continue
    if house3 == 'arnold' and lunches[2] != 'soup':
        continue
    if house4 == 'arnold' and lunches[3] != 'soup':
        continue
    
    # Check if Bob loves stir fry
    if house1 == 'bob' and lunches[0] != 'stir fry':
        continue
    if house2 == 'bob' and lunches[1] != 'stir fry':
        continue
    if house3 == 'bob' and lunches[2] != 'stir fry':
        continue
    if house4 == 'bob' and lunches[3] != 'stir fry':
        continue
    
    # Check if Bob is the Butterscotch smoothie drinker
    if house1 == 'bob' and smoothies[0] != 'butterscotch':
        continue
    if house2 == 'bob' and smoothies[1] != 'butterscotch':
        continue
    if house3 == 'bob' and smoothies[2] != 'butterscotch':
        continue
    if house4 == 'bob' and smoothies[3] != 'butterscotch':
        continue
    
    # Check if the person who loves the soup is directly left of the Darkness smoothie drinker
    if (lunches[0] == 'soup' and smoothies[1] != 'darkness') or \
       (lunches[1] == 'soup' and smoothies[2] != 'darkness') or \
       (lunches[2] == 'soup' and smoothies[3] != 'darkness'):
        continue
    
    # Check if Carol is directly left of the Dragonfruit smoothie lover
    if (house1 == 'carol' and smoothies[1] != 'dragonfruit') or \
       (house2 == 'carol' and smoothies[2] != 'dragonfruit') or \
       (house3 == 'carol' and smoothies[3] != 'dragonfruit'):
        continue
    
    # Check if the person whose child is named Bella is the Dragonfruit smoothie lover
    if (children[0] == 'bella' and smoothies[0] != 'dragonfruit') or \
       (children[1] == 'bella' and smoothies[1] != 'dragonfruit') or \
       (children[2] == 'bella' and smoothies[2] != 'dragonfruit') or \
       (children[3] == 'bella' and smoothies[3] != 'dragonfruit'):
        continue
    
    # Check if the Desert smoothie lover and the Darkness smoothie drinker are next to each other
    if (smoothies[0] == 'desert' and smoothies[1] != 'darkness') and \
       (smoothies[1] == 'desert' and smoothies[2] != 'darkness') and \
       (smoothies[2] == 'desert' and smoothies[3] != 'darkness'):
        continue
    
    # Check if the person who is a pizza lover is directly left of the Dragonfruit smoothie lover
    if (lunches[0] == 'pizza' and smoothies[1] != 'dragonfruit') or \
       (lunches[1] == 'pizza' and smoothies[2] != 'dragonfruit') or \
       (lunches[2] == 'pizza' and smoothies[3] != 'dragonfruit'):
        continue
    
    # If all conditions are satisfied, print the name of the person in house 1
    print(house1)
    break