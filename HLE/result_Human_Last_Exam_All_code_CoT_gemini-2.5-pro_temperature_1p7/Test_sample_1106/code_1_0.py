import math

def calculate_optimized_hemisphere_points():
    """
    Solves for the minimized maximum number of points in a hemisphere
    for n points on a (d-1)-dimensional sphere.
    """
    # n = number of points
    n = 15
    # sphere_dim = dimension of the hypersphere
    sphere_dim = 8
    # d = dimension of the ambient space R^d. For S^k, d = k+1.
    d = sphere_dim + 1

    print(f"Solving for n={n} points on a {sphere_dim}-dimensional sphere (in R^{d}).\n")

    # The problem is to find phi(n,d), the value K that results from an optimal point placement.
    # We established a lower bound K >= ceil(n/2)
    lower_bound = math.ceil(n / 2)
    print(f"A simple analysis shows the number must be at least ceil({n}/2) = {lower_bound}.\n")
    
    print("Using a theorem by D.G. Larman (1971) to find the exact value.")
    
    # Larman's theorem for n being odd:
    # phi(n,d) = (n-1)/2 + gamma((n-1)/2, d-1)
    # where gamma(m,q) = 1 if m <= q.

    if n % 2 != 1:
        print(f"Error: The formula used here is for an odd number of points, but n={n} is even.")
        return

    m = (n - 1) // 2
    q = d - 1

    print(f"The formula is: K = (n-1)/2 + gamma(m, q)")
    print(f"Calculating components:")
    print(f"  m = ({n}-1)/2 = {m}")
    print(f"  q = d-1 = {d}-1 = {q}")
    
    gamma_val = -1
    print("\nEvaluating gamma(m, q):")
    if m <= q:
        gamma_val = 1
        print(f"  Condition m <= q is met ({m} <= {q}).")
        print(f"  Therefore, gamma({m}, {q}) = {gamma_val}.")
    else:
        print(f"  Condition m <= q is not met ({m} > {q}).")
        print("  The simplified version of the theorem does not apply.")
        return

    # Final calculation
    result = m + gamma_val

    print("\nFinal calculation:")
    print("K = m + gamma(m, q)")
    # The final equation as requested:
    print(f"{m} + {gamma_val} = {result}")

calculate_optimized_hemisphere_points()