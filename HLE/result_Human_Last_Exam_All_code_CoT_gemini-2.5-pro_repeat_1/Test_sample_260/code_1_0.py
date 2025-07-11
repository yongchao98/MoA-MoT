import fractions

def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    """
    # Step 1: Define the orders of the cyclic groups corresponding to X1, X2, and X3.
    n1 = 5
    n2 = 8
    n3 = 2
    orders = [n1, n2, n3]
    num_groups = len(orders)

    # Step 2: Explain the method and groups involved.
    print("The space Y is the connected sum of X1, X2, and X3.")
    print(f"The fundamental group is G = pi_1(Y) = Z_{n1} * Z_{n2} * Z_{n3}.")
    print(f"The first homology group is H_1(Y) = Z_{n1} + Z_{n2} + Z_{n3}.")
    print("The kernel K of the Hurewicz map is the commutator subgroup [G, G], which is a free group of unknown rank r.")
    print("We will use the formula relating group Euler characteristics: chi(K) = [G:K] * chi(G).")
    print("-" * 20)

    # Step 3: Perform the calculations.

    # Calculate the index [G:K] = |H_1(Y)|
    index = 1
    for n in orders:
        index *= n

    # Calculate the Euler characteristic of G. We use the fractions module for precision.
    chi_G = fractions.Fraction(0)
    for n in orders:
        chi_G += fractions.Fraction(1, n)
    chi_G -= (num_groups - 1)

    # Calculate the Euler characteristic of K.
    chi_K = index * chi_G

    # The rank r is derived from chi(K) = 1 - r.
    rank = 1 - chi_K

    # Step 4: Print the numerical steps of the calculation.
    print(f"The index [G:K] is the order of H_1(Y):")
    print(f"[G:K] = {n1} * {n2} * {n3} = {index}")
    print("\nThe Euler characteristic of G is calculated as:")
    chi_G_str = " + ".join([f"1/{n}" for n in orders])
    print(f"chi(G) = ({chi_G_str}) - ({num_groups - 1}) = {chi_G.numerator}/{chi_G.denominator}")
    print("\nThe Euler characteristic of the kernel K is:")
    print(f"chi(K) = [G:K] * chi(G) = {index} * ({chi_G.numerator}/{chi_G.denominator}) = {int(chi_K)}")
    print("\nFinally, the rank r of the free group K is found using chi(K) = 1 - r:")
    print("r = 1 - chi(K)")
    print(f"r = 1 - ({int(chi_K)}) = {int(rank)}")

if __name__ == '__main__':
    solve_rank_of_kernel()
    print("<<<95>>>")
