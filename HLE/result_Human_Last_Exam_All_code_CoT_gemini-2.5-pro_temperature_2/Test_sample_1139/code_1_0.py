import math

def calculate_bosonic_variables():
    """
    Calculates the number of non-Grassmann (bosonic) variables for the
    supersymmetric sigma-model in symmetry class D with two replicas.
    """
    # Number of replicas
    n_r = 2
    
    print(f"The number of replicas is n_r = {n_r}.")
    print("The target manifold for the bosonic sector is derived from the symmetric superspace G/H,")
    print("where G = Osp(2*n_r|2*n_r) and H = Osp(n_r|n_r) x Osp(n_r|n_r).")
    print("The number of bosonic variables is dim(G_B) - dim(H_B).\n")

    # Define dimension calculation functions
    def dim_O(k):
        """Calculates the dimension of the Orthogonal group O(k)."""
        return k * (k - 1) // 2

    def dim_Sp(two_k):
        """Calculates the dimension of the Symplectic group Sp(2k)."""
        if two_k % 2 != 0:
            raise ValueError("Argument to Sp must be even.")
        k = two_k // 2
        return k * (2 * k + 1)

    # Calculate dimension of G_B = O(2*n_r) x Sp(2*n_r) = O(4) x Sp(4)
    g_o_k = 2 * n_r
    g_sp_2k = 2 * n_r
    dim_O_G = dim_O(g_o_k)
    dim_Sp_G = dim_Sp(g_sp_2k)
    dim_G_B = dim_O_G + dim_Sp_G
    
    print("1. Calculate the dimension of the group G_B = O(4) x Sp(4):")
    print(f"   dim(O({g_o_k})) = {g_o_k}*({g_o_k}-1)/2 = {dim_O_G}")
    print(f"   dim(Sp({g_sp_2k})) = ({g_sp_2k}/2) * (2*({g_sp_2k}/2) + 1) = {dim_Sp_G}")
    print(f"   dim(G_B) = {dim_O_G} + {dim_Sp_G} = {dim_G_B}\n")

    # Calculate dimension of H_B = [O(n_r) x Sp(n_r)] x [O(n_r) x Sp(n_r)]
    # We calculate for one part and multiply by 2.
    # H_B_part = O(2) x Sp(2)
    h_o_k = n_r
    h_sp_2k = n_r
    dim_O_H = dim_O(h_o_k)
    dim_Sp_H = dim_Sp(h_sp_2k)
    dim_H_B_part = dim_O_H + dim_Sp_H
    dim_H_B = 2 * dim_H_B_part

    print("2. Calculate the dimension of the subgroup H_B = [O(2) x Sp(2)] x [O(2) x Sp(2)]:")
    print(f"   Dimension of one part [O({h_o_k}) x Sp({h_sp_2k})]:")
    print(f"     dim(O({h_o_k})) = {h_o_k}*({h_o_k}-1)/2 = {dim_O_H}")
    print(f"     dim(Sp({h_sp_2k})) = ({h_sp_2k}/2) * (2*({h_sp_2k}/2) + 1) = {dim_Sp_H}")
    print(f"     Sub-dimension = {dim_O_H} + {dim_Sp_H} = {dim_H_B_part}")
    print(f"   Total dim(H_B) = 2 * {dim_H_B_part} = {dim_H_B}\n")

    # Final calculation
    total_vars = dim_G_B - dim_H_B
    
    print("3. Calculate the total number of bosonic variables:")
    print(f"   Number = dim(G_B) - dim(H_B)")
    print(f"   Number = {dim_G_B} - {dim_H_B} = {total_vars}")


if __name__ == '__main__':
    calculate_bosonic_variables()