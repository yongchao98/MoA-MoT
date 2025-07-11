import math

def solve_and_print_sum():
    """
    This function solves the described problem by calculating the sum over all
    natural dimensions of the specified ratio.

    The ratio for each dimension 'd' is R_d = A_d / B_d, where:
    - A_d is the expected volume of a d-simplex with vertices at the origin and
      d points chosen uniformly from orthogonal axes within a unit ball.
      A_d = (1/d!) * (1/2)^d
    - B_d is the expected pairwise Euclidean distance between these d+1 points.
      B_d = (1 + (d-1)*C) / (d+1), where C is a constant.

    The constant C is the expected distance E[sqrt(r_i^2 + r_j^2)] for r_i, r_j
    chosen uniformly from (-1, 1). This is equivalent to the average distance
    between two random points in a unit square (Robbins' constant).
    C = (sqrt(2) + 2 * ln(1 + sqrt(2))) / 3.

    The function calculates the sum S = sum(R_d for d=1 to infinity).
    """
    # Define the constant C (Robbins' constant)
    C = (math.sqrt(2) + 2 * math.log(1 + math.sqrt(2))) / 3

    total_sum = 0.0
    # Summation up to d=30 is sufficient for convergence to the required precision.
    max_dim = 30
    
    equation_terms = []

    for d in range(1, max_dim + 1):
        # The term for dimension d is R_d = (d+1) / (d! * 2^d * (1 + (d-1)*C))
        try:
            term = (d + 1) / (math.factorial(d) * (2**d) * (1 + (d - 1) * C))
        except OverflowError:
            # For large d, the term is effectively zero.
            term = 0.0
        
        total_sum += term
        
        # Store the first 5 terms for the equation string.
        if len(equation_terms) < 5:
            equation_terms.append(f"{term:.3f}")
        
        # Break early if the rest of the series is negligible.
        if term < 1e-15:
            break
    
    # Add ellipsis to indicate an infinite series.
    equation_terms.append("...")

    # Construct the final equation string as requested.
    equation_str = " + ".join(equation_terms)
    
    # Print the final output showing the numbers in the equation.
    print("The final sum is determined by the following equation:")
    print(f"{equation_str} = {total_sum:.3f}")

if __name__ == '__main__':
    solve_and_print_sum()