# Number of eggs in a dozen
dozen_eggs = 12

# Total eggs from dozens and loose eggs
total_eggs = 4 * dozen_eggs + 2

# Eggs needed per glass
eggs_per_glass = 5

# Glasses that can be made
glasses = total_eggs // eggs_per_glass

# Glasses per tray
glasses_per_tray = 5

# Trays that can be put out
trays = glasses // glasses_per_tray

print(trays)