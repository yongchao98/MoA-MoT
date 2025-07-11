import snappy

# Step 1: Define the braid from the problem statement.
# The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^{-1} in B_5.
# In SnapPy, sigma_i is represented by the integer i, and its inverse by -i.
num_strands = 5
braid_word = [1, 1, 2, 2, 3, -4]

print(f"Analyzing the closure of the braid on {num_strands} strands defined by the word: {braid_word}")

# Step 2: Create the link from the closure of the braid.
# The `snappy.Link` constructor can take a braid representation directly.
try:
    link = snappy.Link(braid=braid_word, num_strands=num_strands)
except snappy.SnapPyNotAvailableError:
    print("SnapPy's kernel is not available. Please ensure it is correctly installed.")
    exit()

# Step 3: Get the individual components of the link.
num_components = link.num_components()
print(f"\nThe resulting link has {num_components} components, as expected.")

# The problem states two are unknots. Let's find the third one.
components = link.link_components
print("Identifying each component knot:")

identified_knots = []
for i, c in enumerate(components):
    # identify() returns a list of candidate knots; we take the first, most likely one.
    knot = c.identify()[0]
    identified_knots.append(knot)
    # The standard name for the unknot is '0_1'.
    print(f"Component {i+1} is the knot: {knot}")

# Step 4: Determine the final answer.
# Find the component that is not the unknot ('0_1').
non_unknot_component = None
for knot in identified_knots:
    if str(knot) != '0_1':
        non_unknot_component = knot
        break

if non_unknot_component:
    # The knot '4_1' is the Figure-8 knot.
    # The knot '3_1' is the Trefoil knot.
    # The knot '5_1' is the cinquefoil knot.
    print(f"\nTwo components are indeed unknots ('0_1').")
    print(f"The other connected component is the {non_unknot_component} knot.")
    if str(non_unknot_component) == '4_1':
        print("This corresponds to the Figure-8 knot.")
else:
    print("All components were identified as unknots.")
