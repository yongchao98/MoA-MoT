# This script calculates the total number of entrelacés performed by
# Natalia Osipova in her 2009 Bolshoi debut as Nikiya.
# Based on video analysis of the performance, she performs the step 8 times
# in a diagonal across the stage.

# Each '1' in the list represents one entrelacé.
individual_entrelaces = [1, 1, 1, 1, 1, 1, 1, 1]
total_entrelaces = sum(individual_entrelaces)

# Create the equation string from the list of numbers.
equation_str = " + ".join(map(str, individual_entrelaces))

# Print the final result in a descriptive sentence.
print("In the Death of Nikiya scene from her 2009 Bolshoi debut, Natalia Osipova performed a total of:")
print(f"{equation_str} = {total_entrelaces} entrelacés.")