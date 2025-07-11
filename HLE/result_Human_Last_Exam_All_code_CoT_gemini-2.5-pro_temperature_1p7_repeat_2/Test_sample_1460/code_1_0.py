import spherogram as sp
import sympy
import sys

# Step 0: Explain the task and the braid definition
print("This script analyzes the link formed by the closure of a braid.")
print("The braid group is B_5 (5 strands).")
# As per the prompt's instruction to output the numbers in the equation
print("The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3^1 * sigma_4^-1.")
print("The numbers defining the braid are indices {1, 2, 3, 4} and powers {2, 2, 1, -1}.\n")


# Step 1: Define the braid group and construct the braid
try:
    # Define the braid group on 5 strands
    B5 = sp.BraidGroup(5)
    s = B5.generators()
except ImportError as e:
    print(f"Error: Spherogram library not found. Please install it using 'pip install spherogram'.")
    print(f"Original error: {e}")
    sys.exit(1)


# Braid generators s[i] correspond to sigma_{i+1}
# beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1
beta = s[0]**2 * s[1]**2 * s[2] * s[3]**-1


# Step 2: Take the closure of the braid to form a link
L = beta.closure()


# Step 3: Get the components of the link
components = L.link_components
num_components = len(components)
print(f"Analysis of the link L = closure(beta):")
print(f"The link has {num_components} connected components.\n")


# Step 4: Analyze each component
# Define the variable for the polynomial and the known polynomials for comparison
t = sympy.var('t')
knot_type_map = {
    sympy.Integer(1): "Unknot",
    t - 1 + 1/t: "Trefoil",
    -t + 3 - 1/t: "Figure-8",
    t**2 - t + 1 - t**-1 + t**-2: "5_1 knot"
}

unknot_count = 0
other_knots = []

print("Identifying each component by its Alexander polynomial:")
for i, K in enumerate(components):
    # Compute the (normalized) Alexander polynomial for the component
    poly = K.alexander_polynomial()
    
    knot_name = "Unknown"
    # Identify the knot by its polynomial
    for p, name in knot_type_map.items():
        # The polynomials should match exactly since spherogram normalizes them
        if sympy.simplify(poly - p) == 0:
            knot_name = name
            break
            
    print(f"Component {i+1}:")
    print(f"  - Alexander polynomial: {poly}")
    print(f"  - Identified as: {knot_name}")
    
    if knot_name == "Unknot":
        unknot_count += 1
    else:
        other_knots.append(knot_name)

print("\n-------------------------------------------")
# Step 5: State the final conclusion
print("Conclusion:")
print(f"The analysis found {unknot_count} unknot(s) and {len(other_knots)} other knot(s).")

if unknot_count == 2 and len(other_knots) == 1:
    print("This matches the problem statement that two components are unknots.")
    print(f"The other connected component is the {other_knots[0]} knot.")
elif num_components != 3:
     print("The number of components is not 3, which contradicts the premise.")
else:
    print("The condition that exactly two components are unknots is not met by this braid.")
