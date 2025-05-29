# Total number of animals
total_animals = 3016226

# Total number of legs
total_legs = 70

# Solve the equations
# c + w = total_animals
# c + 2w = total_legs / 2

# Calculate number of cows
w = (total_legs // 2) - total_animals

# Calculate number of chickens
c = total_animals - w

print(c)