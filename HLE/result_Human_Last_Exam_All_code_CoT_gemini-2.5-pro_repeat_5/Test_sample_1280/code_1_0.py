import math

def solve():
    """
    This function explains the reasoning to find the maximum number of differing energy levels
    between two supersymmetric partner Hamiltonians H_0 and H_1.
    """

    # Step 1 & 2: Relationship between H_0 and H_1 spectra
    # H_0 = L^+L - alpha
    # H_1 = LL^+ - alpha
    # The operators L and L^+ map the eigenstates of one Hamiltonian to the other.
    # If H_0 * psi = E * psi, then H_1 * (L*psi) = E * (L*psi).
    # If H_1 * chi = E * chi, then H_0 * (L^+*chi) = E * (L^+*chi).
    # This creates a pairing of energy levels. The spectra can only differ if this mapping fails.

    # Step 3: Source of difference - the zero modes
    # The mapping L*psi fails to produce a valid eigenstate if L*psi = 0.
    # The mapping L^+*chi fails if L^+*chi = 0.
    # These are the "zero modes".
    # If L*psi = 0, then H_0*psi = (L^+L - alpha)*psi = L^+(0) - alpha*psi = -alpha*psi.
    # This unpaired state for H_0 has energy E = -alpha.
    # If L^+*chi = 0, then H_1*chi = (LL^+ - alpha)*chi = L(0) - alpha*chi = -alpha*chi.
    # This unpaired state for H_1 has energy E = -alpha.
    # The total number of differing levels is the sum of the number of zero modes for L and L^+.
    # Number of differing levels = dim(ker(L)) + dim(ker(L^+)).

    # Step 4: Maximum number of zero modes for L and L^+
    # The equation L*psi = (d/dx - W(x))*psi = 0 is a first-order linear ODE.
    # For a given set of physical boundary conditions on [0,1], such an equation
    # can have at most one linearly independent solution.
    # So, dim(ker(L)) can be 0 or 1.
    dim_ker_L_max = 1

    # Similarly, L^+*chi = (-d/dx - W(x))*chi = 0 is also a first-order linear ODE.
    # It can also have at most one linearly independent solution.
    # So, dim(ker(L^+)) can be 0 or 1.
    dim_ker_L_plus_max = 1

    # Step 5: Final calculation
    # The maximum number of differing levels is the maximum possible value of dim(ker(L)) + dim(ker(L^+)).
    max_diff_levels = dim_ker_L_max + dim_ker_L_plus_max

    # This maximum is achievable. For example, with periodic boundary conditions on [0,1]
    # and a periodic superpotential W(x) with integral(W(x)dx) from 0 to 1 = 0,
    # both dim(ker(L)) and dim(ker(L^+)) will be 1.

    # The equation is: Max_Levels = Max_dim(ker(L)) + Max_dim(ker(L^+))
    # We output the numbers in this conceptual equation.
    print(f"The maximum number of levels that can differ is given by the sum of the maximum possible dimensions of the kernels of L and L^+.")
    print(f"Maximum dimension of ker(L) = {dim_ker_L_max}")
    print(f"Maximum dimension of ker(L^+) = {dim_ker_L_plus_max}")
    print(f"The final equation is: {max_diff_levels} = {dim_ker_L_max} + {dim_ker_L_plus_max}")
    print(f"So, the maximum number of differing levels is {max_diff_levels}.")

solve()