import sys

def get_number_of_maximal_tori_e8():
    """
    This function determines the number of Fq-rational maximal tori
    for a reductive group of type E8.

    The solution relies on a well-established theorem in the theory of
    algebraic groups.
    """
    
    # The algebraic group G is of type E8.
    group_type = "E_8"
    
    # The corresponding Weyl group is W(E_8).
    weyl_group_name = "W(E_8)"
    
    # According to a fundamental theorem for reductive groups over finite fields,
    # the number of Fq-rational maximal tori (up to G(Fq)-conjugacy)
    # is equal to the number of conjugacy classes of its Weyl group.
    
    # The number of conjugacy classes for the Weyl group W(E_8) is a
    # known value from the study of reflection groups.
    # The character table of W(E_8), which was famously computed, has 112 rows
    # and 112 columns. The number of columns corresponds to the number of
    # conjugacy classes.
    num_conjugacy_classes = 112
    
    # Therefore, the number of Fq-rational maximal tori is this value.
    num_tori = num_conjugacy_classes
    
    print(f"Let G be a reductive group of type {group_type} over the finite field Fq.")
    print("The number of Fq-rational maximal tori, up to G(Fq)-conjugacy, is equal to the number of conjugacy classes of its Weyl group W.")
    print(f"For type {group_type}, the Weyl group is {weyl_group_name}.")
    
    # Following the prompt to output numbers in the final equation.
    # The "equation" is the equality stated by the theorem.
    print(f"Number of conjugacy classes of {weyl_group_name} = {num_conjugacy_classes}")
    print(f"Final Answer: The number of Fq-rational maximal tori of G is {num_tori}.")

if __name__ == '__main__':
    get_number_of_maximal_tori_e8()