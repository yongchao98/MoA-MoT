# The user needs to have pyknotid installed. It can be installed with:
# pip install pyknotid
import pyknotid.catalogue as cat

# The knot component is equivalent to the closure of a 3-strand braid.
# This braid is represented by the word sigma_1 * sigma_2^(-1).
# In numerical notation for pyknotid, this corresponds to the list of numbers [1, -2].
braid_word_numbers = [1, -2]

# Create the knot object from the braid closure
# from_braid returns a list of components, but our braid has one component
knot = cat.from_braid(braid_word_numbers)[0]

# Identify the knot from the catalogue
# This returns a list of possible matches. We'll take the first one.
identified_knots = knot.identify()
knot_name = "Unknown"
if identified_knots:
    # The name is usually given as a tuple, e.g., (3, 1) for the trefoil knot 3_1
    crossings, index = identified_knots[0]
    if crossings == 3 and index == 1:
        knot_name = "Trefoil"
    elif crossings == 4 and index == 1:
        knot_name = "Figure-8"
    elif crossings == 5 and index == 1:
        knot_name = "$5_1$ knot"
    elif crossings == 0:
        knot_name = "Unknot"
    else:
        knot_name = f"{crossings}_{index} knot"

# Print the final "equation" and the result
print("The final equation to determine the knot is based on its braid representation.")
print(f"The braid consists of generators represented by the numbers: {braid_word_numbers[0]}, {braid_word_numbers[1]}")
print(f"This corresponds to the braid word sigma_1 * sigma_2^(-1).")
print(f"The closure of this braid forms a knot identified as the: {knot_name} knot.")
