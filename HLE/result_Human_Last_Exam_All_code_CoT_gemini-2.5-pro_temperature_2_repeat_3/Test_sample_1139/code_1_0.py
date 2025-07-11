import math

def calculate_bosonic_variables(N):
    """
    Calculates the number of non-Grassmann variables for the supersymmetric sigma-model
    of symmetry class D with N replicas.

    The target space is the supermanifold G/H where:
    G = OSp(4N|4N)
    H = OSp(2N|2N) x OSp(2N|2N)

    The number of bosonic variables is dim(G0) - dim(H0), where G0 and H0 are
    the bosonic subgroups of G and H.
    """

    # --- Dimensions for G = OSp(4N|4N) ---
    # The bosonic subgroup G0 is O(4N) x Sp(4N,R).
    # dim(O(n)) = n(n-1)/2
    # dim(Sp(2k,R)) = k(2k+1)

    # For G0, n_ortho = 4N, k_symp = 2N
    n_ortho_G = 4 * N
    k_symp_G = 2 * N
    
    dim_O_4N = n_ortho_G * (n_ortho_G - 1) / 2
    dim_Sp_4N = k_symp_G * (2 * k_symp_G + 1)
    dim_G0 = dim_O_4N + dim_Sp_4N

    # --- Dimensions for H = OSp(2N|2N) x OSp(2N|2N) ---
    # H0 consists of two identical factors, each being the bosonic subgroup of OSp(2N|2N).
    # For one factor, the bosonic subgroup is O(2N) x Sp(2N,R).
    # Here, n_ortho = 2N, k_symp = N
    n_ortho_H = 2 * N
    k_symp_H = N

    dim_O_2N = n_ortho_H * (n_ortho_H - 1) / 2
    dim_Sp_2N = k_symp_H * (2 * k_symp_H + 1)
    dim_H0_one_factor = dim_O_2N + dim_Sp_2N
    dim_H0 = 2 * dim_H0_one_factor

    # --- Final Calculation ---
    num_variables = dim_G0 - dim_H0

    print("The number of replicas is N = {}.".format(N))
    print("The target supermanifold is G/H = OSp(4N|4N) / (OSp(2N|2N) x OSp(2N|2N)).")
    print("The number of bosonic variables equals the dimension of the bosonic part of G/H.")
    print("\nStep 1: Calculate the dimension of the bosonic subgroup G0 of G = OSp({}|{}).".format(4*N, 4*N))
    print("G0 = O({}) x Sp({}),R)".format(n_ortho_G, n_ortho_G))
    print("dim(G0) = dim(O({})) + dim(Sp({},R))".format(n_ortho_G, n_ortho_G))
    print("dim(G0) = {} + {} = {}".format(int(dim_O_4N), int(dim_Sp_4N), int(dim_G0)))
    
    print("\nStep 2: Calculate the dimension of the bosonic subgroup H0 of H = OSp({}|{}) x OSp({}|{}).".format(2*N, 2*N, 2*N, 2*N))
    print("H0 has two factors. For one factor OSp({}|{}), the bosonic subgroup is O({}) x Sp({},R).".format(2*N, 2*N, n_ortho_H, n_ortho_H))
    print("Dimension of one factor = dim(O({})) + dim(Sp({},R)) = {} + {} = {}".format(n_ortho_H, n_ortho_H, int(dim_O_2N), int(dim_Sp_2N), int(dim_H0_one_factor)))
    print("dim(H0) = 2 * {} = {}".format(int(dim_H0_one_factor), int(dim_H0)))
    
    print("\nStep 3: Calculate the total number of non-Grassmann (bosonic) variables.")
    print("Number of variables = dim(G0) - dim(H0)")
    print("Final Equation: {} - {} = {}".format(int(dim_G0), int(dim_H0), int(num_variables)))

if __name__ == '__main__':
    # Number of replicas as given in the problem
    N = 2
    calculate_bosonic_variables(N)