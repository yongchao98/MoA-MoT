# The order of Old Russian enclitics is determined by a fixed linguistic template.
# This script stores the correct order of the provided enclitics and prints it.

# 1. же (zhe): Emphatic particle
# 2. бы (by): Conditional particle
# 3. мя (mya): Pronominal clitic
# 4. еси (esi): Verbal clitic
# 5. бо (bo): Causal particle
ordered_enclitics = ['же', 'бы', 'мя', 'еси', 'бо']

# Print the final sequence as an "equation" showing the correct order.
# The prompt asked to output each element in the equation, which we do here.
print("The correct order for the enclitics is:")
final_equation = " -> ".join(ordered_enclitics)
print(final_equation)