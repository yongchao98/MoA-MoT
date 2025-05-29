# Let's create a function to solve the puzzle
def solve_puzzle():
    # Initialize the houses with None
    houses = [None] * 4

    # Assign known values based on the clues
    # Bob is the bird keeper and loves brown, so he must be in house 2 or 3
    # Arnold uses Huawei P50, and Alice is directly left of Arnold
    # The person who loves blue is the fish enthusiast and is directly left of the person who uses Google Pixel 6
    # The person who uses Google Pixel 6 loves purple and is the cat lover

    # Let's try to place Bob in house 2
    houses[1] = {'name': 'bob', 'animal': 'bird', 'color': 'brown'}

    # Arnold cannot be in house 1 or 4, so he must be in house 3
    houses[2] = {'name': 'arnold', 'phone': 'huawei p50'}

    # Alice is directly left of Arnold, so she must be in house 2
    houses[0] = {'name': 'alice'}

    # The person who uses Google Pixel 6 loves purple and is the cat lover, and is directly right of the fish enthusiast
    # The fish enthusiast loves blue
    houses[3] = {'phone': 'google pixel 6', 'color': 'purple', 'animal': 'cat'}
    houses[2]['color'] = 'blue'
    houses[2]['animal'] = 'fish'

    # The person who uses Samsung Galaxy S21 is next to Arnold, so they must be in house 1
    houses[0]['phone'] = 'samsung galaxy s21'

    # Fill in the remaining details
    houses[0]['animal'] = 'dog'
    houses[0]['color'] = 'white'

    # Determine the name of the person in house 1
    return houses[0]['name']

# Execute the function and print the result
print(solve_puzzle())