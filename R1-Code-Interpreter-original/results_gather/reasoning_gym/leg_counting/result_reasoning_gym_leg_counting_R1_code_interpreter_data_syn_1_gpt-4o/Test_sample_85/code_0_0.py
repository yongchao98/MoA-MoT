# Number of animals
cows = 12
cats = 2
lions = 12
wasps = 9
shrimps = 10
grasshoppers = 6

# Number of legs per animal
legs_per_cow = 4
legs_per_cat = 4
legs_per_lion = 4
legs_per_wasp = 6
legs_per_shrimp = 10
legs_per_grasshopper = 6

# Total number of legs
total_legs = (cows * legs_per_cow) + (cats * legs_per_cat) + (lions * legs_per_lion) + \
             (wasps * legs_per_wasp) + (shrimps * legs_per_shrimp) + (grasshoppers * legs_per_grasshopper)

print(total_legs)