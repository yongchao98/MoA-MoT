import math

def main():
    """
    Calculates the number of non-Grassman variables for the specified sigma-model.
    """
    # The problem specifies a "two-replica" model, which corresponds to setting the index n=2
    # in the relevant family of symmetric superspaces for class D: OSp(2n|2n)/U(n|n).
    n = 2

    # The bosonic manifold is (O(2n)/U(n)) x (Sp(2n,R)/U(n)).
    # We need to calculate the dimensions of the groups for n=2.
    # The groups are O(2n) -> O(4), U(n) -> U(2), and Sp(2n,R) -> Sp(4,R).

    # Dimension of the Orthogonal group O(m) = m*(m-1)/2
    m_o = 2 * n
    dim_O = m_o * (m_o - 1) // 2

    # Dimension of the Unitary group U(m) = m^2
    m_u = n
    dim_U = m_u**2

    # Dimension of the real Symplectic group Sp(2m,R) = m*(2m+1)
    m_sp = n
    dim_Sp = m_sp * (2 * m_sp + 1)

    # The total number of bosonic variables is the dimension of the manifold:
    # dim(O(2n)/U(n)) + dim(Sp(2n,R)/U(n))
    # = (dim(O(2n)) - dim(U(n))) + (dim(Sp(2n,R)) - dim(U(n)))
    
    dim_part1 = dim_O - dim_U
    dim_part2 = dim_Sp - dim_U
    total_dimension = dim_part1 + dim_part2

    print("The number of non-Grassman variables is the dimension of the bosonic target space.")
    print("For a two-replica (n=2) model of class D, the space is derived from OSp(4|4)/U(2|2).")
    print("Its bosonic part is (O(4)/U(2)) x (Sp(4,R)/U(2)).")
    print("\nCalculation steps:")
    print(f"dim(O(4)) = 4 * (4-1) / 2 = {dim_O}")
    print(f"dim(U(2)) = 2^2 = {dim_U}")
    print(f"dim(Sp(4,R)) = 2 * (2*2 + 1) = {dim_Sp}")
    
    print("\nThe dimension of the total manifold is calculated as:")
    final_equation = f"(dim(O(4)) - dim(U(2))) + (dim(Sp(4,R)) - dim(U(2))) = ({dim_O} - {dim_U}) + ({dim_Sp} - {dim_U}) = {dim_part1} + {dim_part2} = {total_dimension}"
    print(final_equation)

if __name__ == "__main__":
    main()