import math

def solve_high_dimensional_sum():
    """
    Calculates the sum over all natural dimensions of the ratio between the
    expected volume of a random convex hull and the expected pairwise
    Euclidean distance of its points.
    """
    # Step 1: Define the constant C.
    # C is the expected distance between two random points on different orthogonal axes,
    # E[sqrt(r_i^2 + r_j^2)], where r_i, r_j are uniform in (-1, 1).
    # This is calculated as the integral of sqrt(x^2 + y^2) over the square [-1,1]x[-1,1]
    # multiplied by the PDF (1/4).
    # The integral evaluates to (1/3) * (sqrt(2) + ln(1 + sqrt(2))).
    C = (math.sqrt(2) + math.log(1 + math.sqrt(2))) / 3.0

    # Step 2: Initialize variables for the summation.
    # The series converges very quickly due to the n! and 2^n terms in the
    # denominator. Summing up to n=50 is more than sufficient for high precision.
    total_sum = 0.0
    n_factorial = 1.0
    max_n = 50

    print("This script calculates the sum of Ratio(n) for n = 1, 2, ...")
    print("The formula for each term is: Ratio(n) = Numerator(n) / Denominator(n)")
    print(f"Where Numerator(n) = (n+1) and Denominator(n) = (2^n * n! * (1 + (n-1)*C)) with C = {C:.6f}\n")
    print("Individual numbers contributing to the sum:")
    print("-" * 50)

    # Step 3: Loop through dimensions n = 1, 2, 3, ... to compute the sum.
    for n in range(1, max_n + 1):
        if n > 1:
            n_factorial *= n

        # Numerator of the ratio expression for dimension n.
        # This is (n+1) from the formula: Ratio(n) = (n+1) / (2^n * n! * (1 + (n-1)C)).
        num_in_ratio = float(n + 1)
        
        # Denominator of the ratio expression.
        power_of_2 = math.pow(2, n)
        den_in_ratio = power_of_2 * n_factorial * (1 + (n - 1) * C)

        # The term for the current dimension n
        term = num_in_ratio / den_in_ratio
        
        # Add the term to the total sum
        total_sum += term
        
        # As requested, output the numbers in the final equation.
        # We show the first few terms being added to the sum.
        if n <= 5:
            print(f"n={n}: Term = {num_in_ratio} / ({power_of_2:.0f} * {n_factorial:.0f} * (1 + {n-1}*C)) = {term:.8f}")
        elif n == 6:
            print("...")

    print("-" * 50)
    print(f"The converged sum over all dimensions is: {total_sum}")
    print(f"\nThe sum with three-decimal precision is: {total_sum:.3f}")

solve_high_dimensional_sum()
<<<1.250>>>