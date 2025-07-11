import math

def solve_e8_tori_problem():
    """
    Calculates and explains the number of F_q-rational maximal tori
    in a reductive group of type E_8.
    """
    # The degrees of the fundamental invariants for a group of type E_8.
    # These are well-known constants in Lie theory.
    degrees_e8 = [2, 8, 12, 14, 18, 20, 24, 30]

    # The number of conjugacy classes in W(E_8) is a known constant.
    # This is the answer to the user's question.
    num_conjugacy_classes = 112

    print("This program determines the number of F_q-rational maximal tori in a group of type E_8.")
    print("-" * 70)
    
    print("\nStep 1: Theoretical Background\n")
    print("The number of G(F_q)-conjugacy classes of F_q-rational maximal tori in a reductive group G")
    print("over a finite field F_q is given by the number of Frobenius-conjugacy classes in its Weyl group W.")
    print("\nFor a group G of type E_8, the group is always 'split' over any field F_q. This means")
    print("the Frobenius action on the Weyl group W(E_8) is trivial. Therefore, the problem reduces to")
    print("finding the number of ordinary conjugacy classes in the Weyl group W(E_8).")
    
    print("\nStep 2: A Related Calculation: The Order of W(E_8)\n")
    print("While the number of conjugacy classes is a known constant, we can use code to compute a related")
    print("fundamental property of W(E_8): its order. The order of a Weyl group can be calculated by")
    print("taking the product of the degrees of its fundamental invariants.")
    
    print("\nThe degrees of the fundamental invariants for E_8 are:", degrees_e8)
    
    # Calculate the order of W(E_8)
    order_w_e8 = math.prod(degrees_e8)

    # Build the equation string to display the calculation
    equation_str = " * ".join(map(str, degrees_e8))

    print(f"\nThe order of W(E_8) is the product of these degrees:")
    print(f"|W(E_8)| = {equation_str}")
    print(f"|W(E_8)| = {order_w_e8:,}")
    
    print("\nStep 3: The Final Answer\n")
    print("The number of conjugacy classes in W(E_8) is a standard result from the literature on")
    print("reflection groups. This value is constant and does not depend on the field size q.")
    
    print(f"\nThe exact number of F_q-rational maximal tori of G is the number of conjugacy classes of W(E_8), which is {num_conjugacy_classes}.")

if __name__ == "__main__":
    solve_e8_tori_problem()