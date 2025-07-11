import pyknotid.catalogue as cat
import sympy

def solve_knot_problem():
    """
    This function analyzes the closure of the given braid to identify its
    components.
    """
    # 1. Represent the braid
    # The braid is beta = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.
    # In pyknotid, sigma_i is represented by i, and its inverse by -i.
    # The generators are 1-indexed, so sigma_1 is 1, sigma_2 is 2, etc.
    braid_word = [1, 1, 2, 2, 3, -4]
    num_strands = 5
    
    print(f"Analyzing the braid beta = sigma_1^2*sigma_2^2*sigma_3*sigma_4^-1 in B_5.")
    print(f"The braid's numerical representation is: {braid_word} on {num_strands} strands.")
    
    # 2. Compute the link from the braid closure
    print("\nComputing the link components from the braid closure...")
    try:
        link_components = cat.from_braid(braid_word, num_strands=num_strands)
    except ImportError as e:
        print(f"Error: A required dependency for pyknotid might be missing. {e}")
        print("Please ensure numpy, sympy, and networkx are installed ('pip install pyknotid').")
        return
        
    print(f"The closure forms a link with {len(link_components)} components, as expected.")

    # 3. Analyze each component using its Jones polynomial
    print("\nAnalyzing each component:")
    non_unknot_component_poly = None
    
    # Define the variable for our polynomial
    q = sympy.var('q')
    
    for i, component in enumerate(link_components):
        print(f"\n--- Component {i+1} ---")
        
        # Calculate the Jones polynomial
        jp = component.jones_polynomial()
        
        # The problem statement requires outputting the numbers in the final equation.
        # We will print the full polynomial equation for the identified knot.
        print(f"Jones Polynomial: V(q) = {jp}")
        
        if sympy.simplify(jp - 1) == 0:
            print("Result: This component is an Unknot.")
        else:
            print("Result: This component is not an Unknot.")
            non_unknot_component_poly = jp
            
    # 4. Identify the non-unknot component from the answer choices
    print("\n--- Identifying the unknown component ---")
    if non_unknot_component_poly is not None:
        # Known Jones polynomials for comparison
        # Using pyknotid's convention
        trefoil_poly = q + q**3 - q**4
        figure8_poly = q**-2 - q**-1 + 1 - q + q**2
        knot_5_1_poly = 1 - q + q**2 - q**3 + q**4
        
        knot_type = "Unknown"
        final_equation = ""

        if sympy.simplify(non_unknot_component_poly - trefoil_poly) == 0:
            knot_type = "Trefoil"
            final_equation = f"V(q) = {trefoil_poly}"
        elif sympy.simplify(non_unknot_component_poly - figure8_poly) == 0:
            knot_type = "Figure-8"
            final_equation = f"V(q) = {figure8_poly}"
        elif sympy.simplify(non_unknot_component_poly - knot_5_1_poly) == 0:
            knot_type = "5_1 knot"
            final_equation = f"V(q) = {knot_5_1_poly}"
        
        print(f"The non-unknot component's polynomial matches that of the {knot_type} knot.")
        print(f"The final identifying equation is: {final_equation}")
        print("\nThe numbers in this equation are the coefficients and exponents of the polynomial.")
        
    else:
        print("No non-unknot component was found, which contradicts the problem description.")

if __name__ == "__main__":
    solve_knot_problem()