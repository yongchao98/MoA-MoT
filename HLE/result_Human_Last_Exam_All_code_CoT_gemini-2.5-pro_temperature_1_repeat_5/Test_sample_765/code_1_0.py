import math
from scipy.integrate import quad

def solve():
    """
    Calculates the area of the region R analytically by summing the areas
    over the relevant unit squares. This script requires the 'scipy' library.
    You can install it by running 'pip install scipy'.
    """
    
    print("This program calculates the area of the complex region R.")
    print("The solution finds all unit squares [n,n+1)x[m,m+1) that can contain parts of R.")
    print("For each such square, it calculates the area of the region inside and sums them up.")
    
    # Find all unique pairs (n, m) with n<=m where sqrt(n^2+m^2) is an integer.
    # We will use symmetry later to account for pairs where n > m.
    unique_squares_to_calc = []
    all_valid_pairs = []
    for n in range(6):
        for m in range(6):
            val_sqrt = math.sqrt(n**2 + m**2)
            # Check if the square root is an integer
            if val_sqrt == int(val_sqrt):
                k = int(val_sqrt)
                all_valid_pairs.append((n, m))
                if n <= m:
                    # Add to list of squares we need to compute the area for
                    unique_squares_to_calc.append(((n, m), k))
    
    print(f"\nUnit squares (n,m) where sqrt(n^2+m^2) is an integer: {all_valid_pairs}")

    total_area = 0.0
    equation_parts = []
    
    print("\nCalculating area contributions from each unique configuration (where n<=m):")

    # Iterate through the unique squares and calculate the area contribution for each.
    for (n, m), k in unique_squares_to_calc:
        # The area is the integral of the region's height over its width.
        # The region is defined by n<=a<n+1, m<=b<m+1, and a^2+b^2 < (k+1)^2
        R_squared = (k + 1)**2

        def integrand(a):
            # Height of the valid region at a specific 'a'.
            # First, find the upper bound for b from the circle equation a^2+b^2=R^2
            if R_squared < a**2:
                # If 'a' is outside the circle radius, the height is 0.
                return 0
            
            b_upper_from_circle = math.sqrt(R_squared - a**2)
            
            # The height is the length of the segment [m, m+1] that is below b_upper_from_circle.
            # This is max(0, min(b_upper_from_circle, m+1) - m)
            height = max(0, min(b_upper_from_circle, m + 1) - m)
            return height

        # Use scipy.integrate.quad for numerical integration of the height function
        # over the square's width [n, n+1].
        area_contribution, _ = quad(integrand, n, n + 1)
        
        # If n!=m, the pair (m,n) exists and has the same area due to symmetry.
        # So we multiply the contribution by 2.
        if n != m:
            total_area += 2 * area_contribution
            print(f"Area for squares ({n},{m}) and ({m},{n}): 2 * {area_contribution:.4f} = {2 * area_contribution:.4f}")
            equation_parts.append(f"2*{area_contribution:.4f}")
        else: # n == m
            total_area += area_contribution
            print(f"Area for square ({n},{m}): {area_contribution:.4f}")
            equation_parts.append(f"{area_contribution:.4f}")

    print("\nThe total area is the sum of these individual contributions.")
    final_equation_str = " + ".join(equation_parts)
    print(f"Total Area = {final_equation_str}")
    print(f"Final calculated area = {total_area:.4f}")
    print(f"\nThe area of R is approximately {total_area:.2f}")

solve()
<<<11.94>>>