def main():
    """
    Calculates the number of non-Grassman variables needed to parametrize the bosonic sector
    of the supersymmetric sigma-model with two replicas for disordered systems of symmetry class D.
    """

    def dim_O(k):
        """Calculates the dimension of the orthogonal group O(k)."""
        return k * (k - 1) // 2

    def dim_Sp2mR(m):
        """Calculates the dimension of the real symplectic group Sp(2m, R)."""
        return m * (2 * m + 1)

    print("Step 1: Identifying the manifold for symmetry class D with 2 replicas.")
    print("The supersymmetric sigma-model is defined on the superspace G/H = O(4|4) / (O(2|2) x O(2|2)).")
    print("The number of non-Grassman variables is the dimension of the bosonic part of this space,")
    print("which is M_b = [O(4)/(O(2)xO(2))] x [Sp(4,R)/(Sp(2,R)xSp(2,R))].")
    print("-" * 60)

    print("Step 2: Calculating the dimension of the first component space M_O = O(4) / (O(2) x O(2)).")
    k_numerator = 4
    k_denominator = 2
    dim_O_num = dim_O(k_numerator)
    dim_O_den = 2 * dim_O(k_denominator)
    dim_M_O = dim_O_num - dim_O_den
    print(f"Dimension of O({k_numerator}) = {k_numerator}({k_numerator}-1)/2 = {dim_O_num}")
    print(f"Dimension of O({k_denominator}) x O({k_denominator}) = 2 * [{k_denominator}({k_denominator}-1)/2] = {dim_O_den}")
    print(f"Dimension of M_O = {dim_O_num} - {dim_O_den} = {dim_M_O}")
    print("-" * 60)

    print("Step 3: Calculating the dimension of the second component space M_Sp = Sp(4,R) / (Sp(2,R) x Sp(2,R)).")
    m_numerator = 2  # for Sp(2m,R) = Sp(4,R)
    m_denominator = 1  # for Sp(2m,R) = Sp(2,R)
    dim_Sp_num = dim_Sp2mR(m_numerator)
    dim_Sp_den = 2 * dim_Sp2mR(m_denominator)
    dim_M_Sp = dim_Sp_num - dim_Sp_den
    print(f"For Sp(4,R), m={m_numerator}. Dimension = {m_numerator}(2*{m_numerator}+1) = {dim_Sp_num}")
    print(f"For Sp(2,R), m={m_denominator}. Dimension of Sp(2,R) x Sp(2,R) = 2 * [{m_denominator}(2*{m_denominator}+1)] = {dim_Sp_den}")
    print(f"Dimension of M_Sp = {dim_Sp_num} - {dim_Sp_den} = {dim_M_Sp}")
    print("-" * 60)
    
    print("Step 4: Calculating the total number of variables.")
    total_variables = dim_M_O + dim_M_Sp
    print("The total number of non-Grassman variables is the sum of the dimensions of the two component manifolds.")
    print(f"Total Variables = (Dimension of M_O) + (Dimension of M_Sp)")
    print(f"The final equation is: {dim_M_O} + {dim_M_Sp} = {total_variables}")

if __name__ == "__main__":
    main()
