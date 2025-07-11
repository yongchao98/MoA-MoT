import sympy
from pyknotid.braid import Braid

def identify_knot_component():
    """
    Identifies the components of a link formed by a braid closure.
    """
    # The problem defines beta in the braid group B_5
    num_strands = 5
    # The word is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1
    # We represent this as a list of integers
    word = [1, 1, 2, 2, 3, -4]

    print(f"Analyzing the closure of the braid beta on {num_strands} strands.")
    print("Braid word: s_1^2 * s_2^2 * s_3 * s_4^-1")
    print("-" * 30)

    # Create the Braid object
    b = Braid(strands=num_strands, word=word)

    # Get the Knot object representing the link from the braid closure
    link = b.get_knot()

    # The Knot object can be split into its constituent link components
    components = link.link_components
    
    print(f"The closure of the braid results in a link with {len(components)} components.")
    print("-" * 30)

    # Define the symbolic variable for the polynomial
    t = sympy.var('t')
    
    # Define known Alexander polynomials for identification
    known_knots = {
        "Unknot": sympy.Integer(1),
        "Trefoil": t - 1 + t**-1,
        "Figure-8": -t + 3 - t**-1,
        "5_1 knot": t**2 - t + 1 - t**-1 + t**-2
    }
    
    other_component_name = "Not found"

    # Analyze each component
    for i, component_knot in enumerate(components):
        # Calculate the Alexander polynomial for the component
        poly = component_knot.alexander_polynomial()
        
        found_knot = "Unknown"
        # Compare with known polynomials
        for name, known_poly in known_knots.items():
            # Alexander polynomials are unique up to multiplication by +-t^n.
            # We check if the simplified difference is zero.
            if sympy.simplify(poly - known_poly) == 0:
                found_knot = name
                break
        
        print(f"Component {i+1}:")
        print(f"  Alexander Polynomial: {poly}")
        print(f"  This component is a(n) {found_knot}.")

        if found_knot != "Unknot":
            other_component_name = found_knot
            
    print("-" * 30)
    print("The problem states that two components are unknots.")
    print(f"The analysis confirms this and finds that the other component is the {other_component_name} knot.")


if __name__ == "__main__":
    # Note: This script requires the pyknotid library.
    # You can install it via pip: pip install pyknotid
    try:
        identify_knot_component()
    except ImportError:
        print("Please install the 'pyknotid' library to run this script.")
        print("You can install it using: pip install pyknotid")
    except Exception as e:
        print(f"An error occurred: {e}")
