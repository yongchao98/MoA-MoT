# Define the sets based on the premises
writers = set(range(1, 101))  # Let's assume 100 writers
spiders = set(range(1, 201))  # All writers are spiders, and there are additional spiders
lions = set(range(150, 201))  # Some spiders are lions

# Check if there is an intersection between writers and lions
intersection = writers.intersection(lions)

# Print the result
print("Some writers are lions:", len(intersection) > 0)