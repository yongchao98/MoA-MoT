# This script calculates the number of depressions in the lacrimal recess of Irritator challengeri.
# According to the scientific literature (Sues et al., 2002), the lacrimal recess
# is subdivided by a ridge into two smaller fossae (depressions).

# Number of depressions in the dorsal (upper) part of the recess.
dorsal_fossa = 1

# Number of depressions in the ventral (lower) part of the recess.
ventral_fossa = 1

# Calculate the total number of depressions.
total_depressions = dorsal_fossa + ventral_fossa

print("The lacrimal recess in Irritator challengeri contains a dorsal fossa and a ventral fossa.")
print(f"This results in a total number of depressions calculated as follows:")
print(f"{dorsal_fossa} + {ventral_fossa} = {total_depressions}")