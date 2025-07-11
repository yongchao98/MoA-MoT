import numpy as np
from scipy.spatial.distance import directed_hausdorff

def demonstrate_perfectness():
    """
    This function demonstrates that the space P(X) is perfect, meaning it has
    no isolated points. We do this for X = [0, 1].

    We take an example set A from P(X) and construct a sequence of sets B_k,
    such that B_k converges to A, proving A is not isolated.
    """
    print("Demonstrating that P(X) has no isolated points.")
    print("Let X be the unit interval [0, 1].")

    # For computational purposes, we approximate an infinite set from P(X)
    # with a large number of points.
    # Let A = {1/n | n=1,2,3...} U {0}.
    n_points_approximation = 2000
    A_points = {1 / n for n in range(1, n_points_approximation + 1)}
    A_points.add(0.0)

    # Convert set to a NumPy array for distance calculations.
    # The shape must be (n_points, n_dimensions).
    A = np.array(list(A_points)).reshape(-1, 1)

    # We will perturb a point in A, say x_p = 1/10, to create a new set.
    x_p = 0.1

    print(f"\nLet A be the set {{1/n | n=1,...,{n_points_approximation}}} U {{0}}.")
    print(f"We will show A is not isolated by perturbing the point {x_p}.")

    # The sequence of sets B_k is constructed by replacing x_p with a
    # point z_k that gets closer and closer to x_p.
    # Let z_k = x_p - 1/k. For z_k to be different from points in A,
    # we choose k appropriately large.
    k_values = [20, 50, 100, 500, 2000]

    print("-" * 40)
    for k in k_values:
        # z_k approaches x_p as k -> infinity.
        z_k = x_p - 1 / k
        
        # Construct B_k's set of points
        B_k_points = A_points.copy()
        if x_p in B_k_points:
            B_k_points.remove(x_p)
        B_k_points.add(z_k)
        B_k = np.array(list(B_k_points)).reshape(-1, 1)

        # The Hausdorff distance d_H(U, V) is max(d_directed(U,V), d_directed(V,U)).
        # In our construction, d_H(A, B_k) is simply the distance |x_p - z_k| = 1/k.
        # This is because all other points are common to both sets.
        hausdorff_dist = max(directed_hausdorff(A, B_k)[0], directed_hausdorff(B_k, A)[0])

        print(f"For k = {k}, we form B_k by replacing {x_p} with z_k = {z_k:.6f}.")
        print(f"The set B_k is distinct from A and also belongs to P(X).")
        print(f"The Hausdorff distance d_H(A, B_k) = {hausdorff_dist:.6f}, which approaches 0 as k increases.")
        print("")

    print("This shows we can find sets in P(X) arbitrarily close to A.")
    print("Therefore, A is not an isolated point. This logic applies to any point in P(X).")
    
    # Final conclusion based on the mathematical reasoning explained.
    # The question asks to output numbers from the final "equation".
    # The cardinality is C = 2^{\aleph_0}.
    print("\n" + "="*40)
    print("The final answer is the cardinality of the continuum.")
    print("This is expressed as: C = 2 ^ aleph_0")
    print("The numbers in this expression are:")
    print("Base:", 2)
    print("The '0' in aleph_0, indicating the smallest infinite cardinal:", 0)
    print("="*40)

demonstrate_perfectness()
