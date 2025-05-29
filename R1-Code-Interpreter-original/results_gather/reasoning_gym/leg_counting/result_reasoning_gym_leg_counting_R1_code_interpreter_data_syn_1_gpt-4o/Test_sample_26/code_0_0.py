# Number of animals
lions = 6
horses = 5
jellyfish = 13
shrimps = 15
cockroaches = 5
wasps = 8

# Number of legs per animal
legs_per_lion = 4
legs_per_horse = 4
legs_per_jellyfish = 0
legs_per_shrimp = 10
legs_per_cockroach = 6
legs_per_wasp = 6

# Total legs calculation
total_legs = (lions * legs_per_lion +
              horses * legs_per_horse +
              jellyfish * legs_per_jellyfish +
              shrimps * legs_per_shrimp +
              cockroaches * legs_per_cockroach +
              wasps * legs_per_wasp)

print(total_legs)