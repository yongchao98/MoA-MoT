import cypari2
import sys

def calculate_conductor():
    """
    Calculates the conductor of the curve y^2 = f(x) using the cypari2 library.
    The function first defines the polynomial f(x) based on the problem's equation
    and then calls the PARI/GP `hyperellconductor` function to find the conductor.
    """
    print("This script calculates the conductor of a hyperelliptic curve.")
    print("Please ensure the 'cypari2' library is installed (`pip install cypari2`).")
    print("-" * 30)

    try:
        # Initialize the PARI/GP environment
        pari = cypari2.Pari()
    except Exception as e:
        print(f"Error: Could not initialize the cypari2 library.", file=sys.stderr)
        print(f"Details: {e}", file=sys.stderr)
        print("Please make sure PARI/GP is installed on your system and cypari2 is installed correctly.", file=sys.stderr)
        sys.exit(1)

    # The curve is defined by the equation:
    # y^2 = 1*x^6 + 4*x^5 + 6*x^4 + 2*x^3 + 1*x^2 + 2*x + 1

    # The coefficients of the polynomial on the right-hand side are explicitly listed here.
    c6, c5, c4, c3, c2, c1, c0 = 1, 4, 6, 2, 1, 2, 1

    # We construct the polynomial as a string to be passed to PARI/GP.
    poly_str = f"{c6}*x^6 + {c5}*x^5 + {c4}*x^4 + {c3}*x^3 + {c2}*x^2 + {c1}*x + {c0}"
    
    print(f"Curve equation: y^2 = {poly_str}")

    # Convert the Python string into a PARI/GP polynomial object
    f_pari = pari(poly_str)

    # Calculate the conductor using PARI/GP's built-in function `hyperellconductor`.
    print("\nCalculating conductor...")
    conductor = pari.hyperellconductor(f_pari)

    # Print the final result
    print("\nThe conductor of the curve is:")
    print(conductor)

if __name__ == "__main__":
    calculate_conductor()