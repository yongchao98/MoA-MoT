# Define the sets based on the premises
humans = set(['human1', 'human2', 'human3'])
dolphins = set(['human1', 'human2', 'human3', 'dolphin1', 'dolphin2'])
dogs = set(['dolphin1'])  # Only dolphin1 is a dog, based on the premise

# Check if some humans are not dogs
some_humans_not_dogs = any(human not in dogs for human in humans)

print(some_humans_not_dogs)