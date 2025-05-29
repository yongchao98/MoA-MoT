# Number of animals
deers = 9
beetles = 2
giraffes = 15
fireflies = 3
ants = 6
chickens = 14
woodlouses = 2
cats = 6
humans = 6
scorpions = 4

# Number of legs per animal
legs_per_deer = 4
legs_per_beetle = 6
legs_per_giraffe = 4
legs_per_firefly = 6
legs_per_ant = 6
legs_per_chicken = 2
legs_per_woodlouse = 14
legs_per_cat = 4
legs_per_human = 2
legs_per_scorpion = 8

# Total number of legs
total_legs = (deers * legs_per_deer +
              beetles * legs_per_beetle +
              giraffes * legs_per_giraffe +
              fireflies * legs_per_firefly +
              ants * legs_per_ant +
              chickens * legs_per_chicken +
              woodlouses * legs_per_woodlouse +
              cats * legs_per_cat +
              humans * legs_per_human +
              scorpions * legs_per_scorpion)

print(total_legs)