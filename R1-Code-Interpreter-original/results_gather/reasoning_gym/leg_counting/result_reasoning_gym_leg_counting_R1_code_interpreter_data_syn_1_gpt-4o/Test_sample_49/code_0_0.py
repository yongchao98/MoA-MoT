# Number of animals
num_ants = 13
num_shrimps = 14
num_lobsters = 10
num_bees = 1
num_fireflies = 15
num_sea_slugs = 6
num_chickens = 5
num_leeches = 12
num_elephants = 6
num_woodlice = 4

# Legs per animal
legs_per_ant = 6
legs_per_shrimp = 10
legs_per_lobster = 10
legs_per_bee = 6
legs_per_firefly = 6
legs_per_sea_slug = 0
legs_per_chicken = 2
legs_per_leech = 0
legs_per_elephant = 4
legs_per_woodlouse = 14

# Total legs calculation
total_legs = (
    num_ants * legs_per_ant +
    num_shrimps * legs_per_shrimp +
    num_lobsters * legs_per_lobster +
    num_bees * legs_per_bee +
    num_fireflies * legs_per_firefly +
    num_sea_slugs * legs_per_sea_slug +
    num_chickens * legs_per_chicken +
    num_leeches * legs_per_leech +
    num_elephants * legs_per_elephant +
    num_woodlice * legs_per_woodlouse
)

print(total_legs)