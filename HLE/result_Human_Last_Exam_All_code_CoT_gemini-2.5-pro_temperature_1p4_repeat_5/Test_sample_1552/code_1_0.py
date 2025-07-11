def find_second_heat_kernel_coefficient():
    """
    Calculates the numerical coefficient of the Tr(F^2) term in the a_2
    Seeley-DeWitt coefficient for a massless gauged Dirac spinor field.
    """

    # A 4D massless Dirac spinor can be decomposed into two Weyl spinors (left-handed and right-handed).
    # The calculation is performed by summing the contributions from each Weyl spinor.
    num_weyl_spinors = 2

    # The coefficient of the Tr(F^2) term in a_2 for a single Weyl spinor is a
    # well-established result from quantum field theory.
    weyl_coeff = -1/6

    # The total coefficient for the Dirac spinor is the sum of the contributions
    # from its constituent Weyl spinors.
    dirac_coeff = num_weyl_spinors * weyl_coeff

    print("This script calculates the numerical coefficient 'C' of the Yang-Mills term")
    print("C * ∫ Tr(F_{\mu\nu}F^{\mu\nu}) dvol")
    print("within the second heat kernel coefficient (a_2) for a massless gauged Dirac spinor.")
    print("\nThe approach is to sum the contributions from the two Weyl spinors that compose a Dirac spinor.")
    print("\nThe final equation is:")
    print(f"{num_weyl_spinors} * ({weyl_coeff}) = {dirac_coeff}")

if __name__ == '__main__':
    find_second_heat_kernel_coefficient()