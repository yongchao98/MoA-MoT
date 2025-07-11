import math

def comb(n, k):
    """
    Computes the binomial coefficient "n choose k".
    Returns 0 if k < 0 or k > n.
    """
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def get_mu_j(j, v=50, k=5):
    """
    Calculates the eigenvalue mu_j of the adjacency matrix A_2 for the Johnson scheme J(v,k).
    The formula is from Brouwer, Cohen, and Neumaier, "Distance-Regular Graphs" (1989), p. 255.
    This corresponds to the graph where vertices are adjacent if their intersection size is k-2.
    """
    i = 2
    val = 0
    # The eigenvalue of A_i is mu_i(j)
    # mu_i(j) = sum_{l=0 to i} (-1)^l * comb(j,l) * comb(k-j, i-l) * comb(v-k-j+l, i-l)
    term1 = comb(j,0) * comb(k-j, 2) * comb(v-k-j, 2)
    term2 = comb(j,1) * comb(k-j, 1) * comb(v-k-j+1, 1)
    term3 = comb(j,2) * comb(k-j, 0) * comb(v-k-j+2, 0)
    
    return term1 - term2 + term3

def main():
    """
    Main function to execute the plan and find the max order of the Smith group.
    """
    v = 50
    k = 5
    
    # Calculate eigenvalues of A_2
    mu = [get_mu_j(j, v, k) for j in range(k + 1)]
    
    # Number of vertices
    n = comb(v, k)
    
    # Calculate eigenvalues of the Seidel matrix S
    theta = []
    # Principal eigenvalue
    theta_0 = n - 1 - 2 * mu[0]
    theta.append(theta_0)
    # Other eigenvalues
    for j in range(1, k + 1):
        theta_j = -1 - 2 * mu[j]
        theta.append(theta_j)

    abs_theta = [abs(t) for t in theta]

    print("The distinct eigenvalues of the Seidel matrix S are:")
    for i in range(len(theta)):
        print(f"theta_{i} = {theta[i]}")
    print("\nTheir absolute values are:")
    for i in range(len(abs_theta)):
        print(f"|theta_{i}| = {abs_theta[i]}")
        
    # Find the gcd of the absolute eigenvalues
    if len(abs_theta) < 2:
        g = abs_theta[0]
    else:
        g = math.gcd(abs_theta[0], abs_theta[1])
        for i in range(2, len(abs_theta)):
            g = math.gcd(g, abs_theta[i])

    print(f"\nThe greatest common divisor (g) of these values is: {g}")

    # The maximum order is the lcm of |theta_j|/g.
    # If g=1, this is just the lcm of |theta_j|.
    
    if g != 1:
        lcm_args = [val // g for val in abs_theta]
    else:
        lcm_args = abs_theta

    if len(lcm_args) < 2:
        max_order = lcm_args[0]
    else:
        max_order = math.lcm(lcm_args[0], lcm_args[1])
        for i in range(2, len(lcm_args)):
            max_order = math.lcm(max_order, lcm_args[i])

    print("\nThe maximum order is the least common multiple of these absolute values (divided by g if g > 1).")
    lcm_str = ", ".join(map(str, lcm_args))
    print(f"Max order = lcm({lcm_str})")
    print(f"Result = {max_order}")
    print("\nThe final equation for the maximum order is:")
    final_eq = f"lcm({', '.join(map(str, abs_theta))})"
    print(f"{final_eq} = {max_order}")


if __name__ == "__main__":
    main()