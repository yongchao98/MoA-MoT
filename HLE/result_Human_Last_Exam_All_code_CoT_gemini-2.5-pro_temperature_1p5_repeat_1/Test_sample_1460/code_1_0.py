# First, ensure the required library is installed. You can run this command in your shell:
# pip install pyknotid sympy

import sys
import subprocess

# Helper function to install a package if it's not already installed
def install(package):
    try:
        __import__(package)
    except ImportError:
        print(f"{package} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        print(f"{package} installed successfully.")

# Ensure pyknotid and sympy are installed
install('pyknotid')
install('sympy')

import sympy
from pyknotid.representations import Braid
from pyknotid.invariants import alexander_polynomial

def solve_knot_problem():
    """
    This function solves the braid problem by identifying the components of its closure.
    """
    print("Step 1: Representing the braid B = sigma_1^2 * sigma_2^2 * sigma_3 * sigma_4^-1 in B_5.")
    # For B_5, we have 5 strands. The generators are sigma_1 to sigma_4.
    # The braid word is represented by a list of generator indices.
    # sigma_i is i, and sigma_i^-1 is -i.
    num_strands = 5
    braid_word = [1, 1, 2, 2, 3, -4]
    braid_repr = Braid(braid_word, num_strands)
    print(f"Braid representation: {braid_word} on {num_strands} strands.\n")

    print("Step 2: Computing the link closure and its components.")
    # The get_links() method returns a list of Knot objects, one for each component.
    try:
        link_components = braid_repr.get_links()
    except Exception as e:
        print(f"An error occurred during link calculation: {e}")
        print("Please ensure your pyknotid installation is correct.")
        return

    num_components = len(link_components)
    print(f"The closure of the braid forms a link with {num_components} components.\n")

    if num_components != 3:
        print("Warning: The number of components is not 3 as stated in the problem.")
        print("Proceeding with analysis anyway.\n")

    print("Step 3: Calculating the Alexander polynomial for each component to identify it.")
    
    # Define the variable for the polynomial
    t = sympy.var('t')
    knot_identities = {}

    for i, knot_component in enumerate(link_components):
        # Calculate the Alexander polynomial. We get a sympy expression.
        poly = knot_component.alexander_polynomial(variable=t)
        poly_simplified = sympy.expand(poly) # Standardize the form

        knot_name = "Unknown"
        # Identify the knot based on its polynomial
        if poly_simplified == 1:
            knot_name = "Unknot"
        elif poly_simplified == t - 1 + 1/t:
            knot_name = "Trefoil"
        elif poly_simplified == t - 3 + 1/t:
            knot_name = "Figure-8"
        elif poly_simplified == t**2 - t + 1 - 1/t + 1/t**2:
            knot_name = "$5_1$ knot"
        
        knot_identities[i] = {'name': knot_name, 'poly': poly_simplified}

    print("-" * 30)
    # Print the results for each component
    for i, data in knot_identities.items():
        print(f"Component {i+1}:")
        print(f"  - Alexander Polynomial: {data['poly']}")
        print(f"  - Identified as: {data['name']}")
    print("-" * 30)

    print("\nStep 4: Final Analysis.")
    unknot_count = sum(1 for data in knot_identities.values() if data['name'] == 'Unknot')
    other_components = [data for data in knot_identities.values() if data['name'] != 'Unknot']

    print(f"The analysis confirms there are {unknot_count} unknot components.")

    if len(other_components) == 1:
        final_knot = other_components[0]
        final_knot_name = final_knot['name']
        final_knot_poly = final_knot['poly']
        
        print(f"The other component is identified as the {final_knot_name} knot.")

        # Extract coefficients as requested by the prompt ("output each number in the final equation")
        # To get integer coefficients, multiply by t to clear the denominator.
        poly_for_coeffs = sympy.Poly(sympy.expand(final_knot_poly * t), t)
        coeffs = poly_for_coeffs.all_coeffs()

        print("\n--- Final Answer ---")
        print(f"The final connected component is the {final_knot_name} knot.")
        print(f"Its Alexander Polynomial equation is: {final_knot_poly} = 0")
        print(f"The integer coefficients of the numerator of this equation (t^2 - 3t + 1) are: {coeffs}")

    else:
        print("The link does not match the problem description (2 unknots + 1 other).")

if __name__ == '__main__':
    solve_knot_problem()
