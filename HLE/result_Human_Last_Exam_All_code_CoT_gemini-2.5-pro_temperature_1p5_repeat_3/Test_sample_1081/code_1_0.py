import math

def calculate_e8_tori_number():
    """
    Calculates the number of F_q-rational maximal tori for a group of type E8.

    This number is equal to the order of the Weyl group of type E_8, |W(E_8)|.
    The order is calculated as the product of the degrees of the fundamental
    invariants of the Weyl group.
    """
    
    # The degrees of the fundamental invariants for the Weyl group of type E_8.
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]

    # Calculate the order of the Weyl group by taking the product of the degrees.
    weyl_group_order = math.prod(degrees)

    # Explanation for the user
    print("The number of F_q-rational maximal tori of a reductive group G of type E_8 is given by the order of its Weyl group, |W(E_8)|.")
    print("A group of type E_8 is simple and simply-connected, so this theorem applies regardless of the field q.")
    print("\nThe order of the Weyl group is the product of the degrees of its fundamental invariants.")
    print("For E_8, the degrees are: 2, 8, 12, 14, 18, 20, 24, 30.")
    
    # Create a string representation of the multiplication
    equation_str = " * ".join(map(str, degrees))
    
    # Print the final calculation and result
    print("\nThe total number of F_q-rational maximal tori is:")
    print(f"|W(E_8)| = {equation_str} = {weyl_group_order:,}")

if __name__ == "__main__":
    calculate_e8_tori_number()