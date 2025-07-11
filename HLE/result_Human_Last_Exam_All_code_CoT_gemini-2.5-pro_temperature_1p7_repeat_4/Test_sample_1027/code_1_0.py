import math

def main():
    """
    This script calculates the dimension of the homology of the given group G
    in degree 31. The problem reduces to calculating the dimension of the
    homology of the abelian group Z^2, which is given by a binomial coefficient.
    """

    # The degree of the homology group
    n = 31
    # The rank of the abelian group Z^k to which the homology of G is isomorphic.
    # A is isomorphic to Z^2, so k=2.
    k = 2

    # The dimension of H_n(Z^k) is given by the binomial coefficient C(k, n).
    # We use math.comb for this calculation.
    # For n > k, math.comb(k, n) correctly returns 0.
    dimension = math.comb(k, n)
    
    # As per the instructions, output the numbers in the final equation.
    print(f"The dimension of H_{n}(G, R) is equal to the dimension of H_{n}(Z^{k}, R).")
    print(f"dim H_{n}(Z^{k}, R) = C(k, n)")
    print(f"dim H_{31}(Z^2, R) = C({k}, {n}) = {dimension}")

if __name__ == "__main__":
    main()
