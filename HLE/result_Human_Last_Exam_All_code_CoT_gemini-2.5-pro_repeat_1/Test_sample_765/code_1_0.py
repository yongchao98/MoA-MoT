import math
from scipy.integrate import quad

def calculate_area_contribution(na, nb):
    """
    Calculates the area of the valid region for a given unit square [na, na+1) x [nb, nb+1).
    """
    m_sq = na*na + nb*nb
    m = math.sqrt(m_sq)
    
    # The condition for non-zero area is that na^2 + nb^2 is a perfect square.
    if m != int(m):
        return 0.0
    m = int(m)
    
    m_plus_1_sq = (m + 1)**2

    def integrand(a):
        """
        This function represents the length of the valid segment in the b-direction for a given a.
        """
        # Calculate bounds for b^2 from the annulus definition m^2 <= a^2 + b^2 < (m+1)^2
        b_sq_low = m_sq - a*a
        b_sq_high = m_plus_1_sq - a*a

        # If b_sq_high is negative, there's no valid b for this a.
        if b_sq_high <= 0:
            return 0.0
        
        # Convert b^2 bounds to b bounds (b must be non-negative)
        b_low_annulus = math.sqrt(b_sq_low) if b_sq_low > 0 else 0.0
        b_high_annulus = math.sqrt(b_sq_high)
        
        # The final range for b is the intersection of the square's side [nb, nb+1]
        # and the annulus's radial range [b_low_annulus, b_high_annulus].
        lower_b = max(nb, b_low_annulus)
        upper_b = min(nb + 1, b_high_annulus)
        
        # The length of the valid segment for b is (upper_b - lower_b), if positive.
        length = upper_b - lower_b
        return max(0.0, length)

    # Integrate the length function over 'a' from na to na+1 to get the area
    area, _ = quad(integrand, na, na + 1)
    return area

def main():
    """
    Main function to find all valid regions, calculate their areas, and sum them up.
    """
    # List of unique pairs (na, nb) with na <= nb that yield valid regions
    # within the 6x6 square (na, nb from 0 to 5).
    pairs_to_check = [
        (0, 0),
        (0, 1), (0, 2), (0, 3), (0, 4), (0, 5),
        (3, 4)
    ]
    
    total_area = 0
    area_terms = []
    
    print("Calculating area contributions for each valid region:")
    
    for na, nb in pairs_to_check:
        area = calculate_area_contribution(na, nb)
        
        if na == nb:
            # Symmetric pair like (0,0) has only one corresponding square
            term = area
            print(f"Area for square ({na},{nb}): {term:.4f}")
        else:
            # Asymmetric pair like (0,1) corresponds to two squares, (0,1) and (1,0)
            term = 2 * area
            print(f"Area for squares ({na},{nb}) and ({nb},{na}): {term:.4f}")
            
        total_area += term
        area_terms.append(term)
        
    print("\nFinal equation with each calculated area term:")
    # Format the equation string. Using .4f for intermediate precision.
    equation_str = " + ".join([f"{term:.4f}" for term in area_terms])
    # The final result is requested to two decimals.
    print(f"{equation_str} = {total_area:.2f}")

    print(f"\nThe total area of R is: {total_area:.2f}")
    
    # Final answer in the required format
    print(f"<<<{total_area:.2f}>>>")

if __name__ == "__main__":
    try:
        main()
    except ImportError:
        print("This script requires the 'scipy' library.")
        print("Please install it using: pip install scipy")
        print("<<<Error: scipy not found>>>")
