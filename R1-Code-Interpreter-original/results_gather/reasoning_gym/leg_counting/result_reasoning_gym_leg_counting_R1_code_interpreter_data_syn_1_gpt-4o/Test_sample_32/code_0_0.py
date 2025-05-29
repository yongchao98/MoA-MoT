# Number of animals
num_humans = 13
num_wasps = 13
num_lobsters = 4
num_sea_slugs = 12
num_lions = 14
num_deers = 13
num_butterflies = 11
num_praying_mantises = 3
num_shrimps = 2
num_insects = 14

# Number of legs per animal
legs_human = 2
legs_wasp = 6
legs_lobster = 10
legs_sea_slug = 0
legs_lion = 4
legs_deer = 4
legs_butterfly = 6
legs_praying_mantis = 6
legs_shrimp = 10
legs_insect = 6

# Total legs calculation
total_legs = (num_humans * legs_human +
              num_wasps * legs_wasp +
              num_lobsters * legs_lobster +
              num_sea_slugs * legs_sea_slug +
              num_lions * legs_lion +
              num_deers * legs_deer +
              num_butterflies * legs_butterfly +
              num_praying_mantises * legs_praying_mantis +
              num_shrimps * legs_shrimp +
              num_insects * legs_insect)

print(total_legs)