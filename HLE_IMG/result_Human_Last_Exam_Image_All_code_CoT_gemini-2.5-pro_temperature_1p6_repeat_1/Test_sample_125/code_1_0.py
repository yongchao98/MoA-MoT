# The plan is to find the maximum numerical locant in the reactant names
# and subtract it from the maximum numerical locant in the product name.

# Reactant locants were extracted from '1,4-difluoro-2-methylbenzene' and '2-acetylnaphthalene'.
reactant_locants = [1, 4, 2, 2]
max_reactant_locant = max(reactant_locants)

# Product locants were extracted from 'as-indaceno[3,2,1,8,7,6-pqrstuv]picene'.
product_locants = [3, 2, 1, 8, 7, 6]
max_product_locant = max(product_locants)

# The number of steps is the difference.
steps = max_product_locant - max_reactant_locant

# Print the final equation with all its numbers.
print(f"{max_product_locant} - {max_reactant_locant} = {steps}")