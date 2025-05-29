# Define sets
musicians = set(range(1, 11))  # Let's assume 10 musicians
whales = set(range(1, 21))     # Let's assume 20 whales
birds = set(range(11, 21))     # Let's assume 10 birds, none of which are musicians

# All musicians are whales
musicians_as_whales = musicians.intersection(whales)

# Some whales are not birds
whales_not_birds = whales.difference(birds)

# Check if some musicians are not birds
musicians_not_birds = musicians_as_whales.intersection(whales_not_birds)

# Output the result
print("Some musicians are not birds:", len(musicians_not_birds) > 0)