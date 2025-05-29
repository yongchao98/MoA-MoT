# Number of animals
elephants = 10
starfishs = 11
cow = 1
butterflys = 2
dogs = 12
shrimps = 11

# Number of legs per animal
legs_per_elephant = 4
legs_per_starfish = 5
legs_per_cow = 4
legs_per_butterfly = 6
legs_per_dog = 4
legs_per_shrimp = 10

# Total number of legs
total_legs = (elephants * legs_per_elephant +
              starfishs * legs_per_starfish +
              cow * legs_per_cow +
              butterflys * legs_per_butterfly +
              dogs * legs_per_dog +
              shrimps * legs_per_shrimp)

print(total_legs)