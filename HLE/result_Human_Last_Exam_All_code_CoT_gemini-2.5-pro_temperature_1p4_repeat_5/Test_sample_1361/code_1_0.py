# The famous sequence of jumps in the Pizzicato solo is set to a standard 8-count musical phrase.
total_counts_in_phrase = 8

# This 8-count phrase is divided evenly between two different types of jumping steps.
# The first half is gargouillades, and the second half is another step (e.g., jetés).
number_of_step_types = 2

# To find the number of gargouillades, we divide the total counts by the number of step types.
number_of_gargouillades = total_counts_in_phrase // number_of_step_types

print(f"In the Act III 'Pizzicato' solo, a key sequence is performed over an 8-count phrase.")
print(f"The phrase is divided into 2 different sets of steps.")
print(f"The number of gargouillades performed is calculated as:")
print(f"{total_counts_in_phrase} / {number_of_step_types} = {number_of_gargouillades}")