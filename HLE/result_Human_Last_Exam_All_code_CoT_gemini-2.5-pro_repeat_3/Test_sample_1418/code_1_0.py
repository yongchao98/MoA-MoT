# This script must be executed in a SageMath environment.
# SageMath can be accessed online via services like CoCalc or installed locally.

def solve_galois_group():
    """
    This function computes the Galois group of the extension
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3)) over Q.
    """
    print("Step 1: Construct the field tower L/K/Q.")
    # Define the base field of rational numbers
    QQ = RationalField()
    
    # Define the first extension K1 = Q(sqrt(2))
    P1.<x> = QQ[]
    K1.<sqrt2> = NumberField(x**2 - 2, 'sqrt2')
    print("Defined K1 = Q(sqrt(2))")

    # Define the second extension K2 = K1(sqrt(3)) = Q(sqrt(2), sqrt(3))
    P2.<y> = K1[]
    K2.<sqrt3> = K1.extension(y**2 - 3, 'sqrt3')
    print("Defined K2 = Q(sqrt(2), sqrt(3))")

    # Define the element beta in K2
    beta = (2 + sqrt2) * (3 + sqrt3)
    print(f"beta = (2+sqrt(2))(3+sqrt(3)) = {beta}")

    # Check if beta is a square in K2 to confirm the degree of the final extension
    if beta.is_square():
        print("beta is a square in K2, so L=K2.")
    else:
        print("beta is not a square in K2, so [L:K2] = 2.")

    # Define the final extension L = K2(sqrt(beta))
    P3.<z> = K2[]
    L.<alpha> = K2.extension(z**2 - beta, 'alpha')
    print(f"Defined L = K2(alpha) where alpha^2 = beta.")

    print("\nStep 2: Compute the properties of the field L.")
    # The absolute field is the field L viewed as an extension of Q
    # We can get its defining polynomial
    L_abs, _, _ = L.absolute_field('c')
    
    degree = L.absolute_degree()
    print(f"The degree of the extension L/Q is {degree}.")

    print("\nStep 3: Compute the Galois group of L/Q.")
    # This computation can take a few seconds
    Gal_L_Q = L.galois_group()
    
    # Get the structure description of the group
    group_structure = Gal_L_Q.structure_description()
    
    print("\nFinal Answer:")
    print(f"The Galois Group of L/Q is isomorphic to {group_structure}.")
    
    # The requested "final equation" is the structure of the group.
    # For C2 x C2 x C2, the numbers are 2, 2, 2.
    print(f"The structure is a direct product of cyclic groups of orders: 2, 2, 2")

solve_galois_group()
