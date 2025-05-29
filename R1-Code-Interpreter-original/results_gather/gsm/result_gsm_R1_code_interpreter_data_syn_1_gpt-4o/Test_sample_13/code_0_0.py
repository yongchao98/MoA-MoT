# Number of eggs in dozens and loose
dozen_eggs = 4
loose_eggs = 2

# Total eggs
total_eggs = dozen_eggs * 12 + loose_eggs

# Eggs per glass and glasses per tray
eggs_per_glass = 5
glasses_per_tray = 5

# Calculate number of glasses and trays
total_glasses = total_eggs // eggs_per_glass
total_trays = total_glasses // glasses_per_tray

print(total_trays)