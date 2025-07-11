import sympy
from pyknotid.make import Braid
from pyknotid.catalogue.manager import Knot

def get_conway_coeff(poly, term):
    """
    Safely extracts the coefficient of a given term from a sympy polynomial.
    Handles both sympy.Poly and general sympy expressions.
    """
    # The variable of the polynomial, typically z for Conway
    z = term.free_symbols.pop()
    
    # Ensure poly is a sympy.Poly object for robust coefficient extraction
    if not isinstance(poly, sympy.Poly):
        poly = sympy.Poly(poly, z)
        
    coeff = poly.coeff_monomial(term)
    return int(coeff)

def main():
    """
    Main function to perform the knot polynomial calculation.
    """
    try:
        # Define a symbolic variable 'z' for our polynomials
        z = sympy.var('z')

        # Step 1: Represent the braid and create the corresponding knot.
        print("--- Knot 1: Closure of Braid \u03B2 ---")
        num_strands = 5
        braid_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
        b = Braid(braid_word, num_strands=num_strands)
        knot_from_braid = b.closing()
        print(f"Braid \u03B2 is represented by: {braid_word} on {num_strands} strands.")

        # Step 2: Calculate the Conway polynomial for the closure of the braid.
        conway_poly_braid = knot_from_braid.conway_polynomial()
        print(f"The Alexander-Conway polynomial for the closure of \u03B2 is: \u2207(z) = {sympy.pretty(conway_poly_braid)}")
        
        # Step 3: Extract the z^2 coefficient from the first polynomial.
        coeff_braid = get_conway_coeff(conway_poly_braid, z**2)
        print(f"The coefficient of z^2 for \u2207(\u03B2) is: {coeff_braid}")

        print("\n--- Knot 2: 10_4 ---")
        # Step 4: Create the knot 10_4 from the Rolfsen table.
        knot_10_4 = Knot.from_table(10, 4)
        print("Fetching knot 10_4 from the standard Rolfsen table.")
        
        # Step 5: Calculate its Conway polynomial.
        conway_poly_10_4 = knot_10_4.conway_polynomial()
        print(f"The Alexander-Conway polynomial for 10_4 is: \u2207(z) = {sympy.pretty(conway_poly_10_4)}")

        # Step 6: Extract the z^2 coefficient.
        coeff_10_4 = get_conway_coeff(conway_poly_10_4, z**2)
        print(f"The coefficient of z^2 for \u2207(10_4) is: {coeff_10_4}")

        # Step 7: Calculate and print the final difference.
        print("\n--- Final Calculation ---")
        difference = coeff_braid - coeff_10_4
        print("The difference in the z^2 coefficients is calculated as:")
        print(f"Coefficient(\u2207(\u03B2)) - Coefficient(\u2207(10_4)) = {coeff_braid} - {coeff_10_4} = {difference}")

    except ImportError:
        print("Please install required libraries: pip install pyknotid sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
<<<0>>>