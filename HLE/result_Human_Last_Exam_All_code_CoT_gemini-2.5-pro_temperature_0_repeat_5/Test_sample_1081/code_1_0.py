def solve_e8_tori_problem():
    """
    Calculates the number of F_q-rational maximal tori of a reductive group of type E_8.

    The solution is based on a known theorem in algebraic group theory.
    """

    # The number of conjugacy classes for the Weyl group W(E_8).
    # This is a known mathematical constant.
    num_conjugacy_classes_W_E8 = 112

    print("Let G be a reductive group of type E_8 over the finite field F_q.")
    print("According to a theorem by Steinberg, the number of F_q-rational maximal tori of G")
    print("is equal to the number of conjugacy classes of its Weyl group, W(E_8).")
    print("\nThis number is a known constant from the theory of reflection groups.")
    
    print("\nThe final equation is:")
    print(f"Number of F_q-rational maximal tori = Number of conjugacy classes of W(E_8) = {num_conjugacy_classes_W_E8}")

if __name__ == "__main__":
    solve_e8_tori_problem()