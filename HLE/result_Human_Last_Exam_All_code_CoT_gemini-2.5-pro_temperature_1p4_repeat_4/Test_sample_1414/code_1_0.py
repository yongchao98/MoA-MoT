import sys
import subprocess

# This script requires the 'mpmath' library.
# The following code tries to install it if it's not already installed.
try:
    import mpmath
except ImportError:
    print("The 'mpmath' library is not installed. Attempting to install...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "mpmath"])
        import mpmath
        print("Successfully installed 'mpmath'.")
    except Exception as e:
        print(f"Failed to install 'mpmath'. Please install it manually using 'pip install mpmath'.")
        print(f"Error: {e}")
        sys.exit(1)

def solve_machin_like_formula():
    """
    Finds the integer coefficients for the given Machin-like formula for pi.
    """
    # Set the precision for the calculations. 100 digits should be sufficient.
    mpmath.mp.dps = 100

    # These are the denominators of the fractions inside the arctan functions.
    denominators = [122, 239, 682, 1252, 2855, 12943]

    # Create the vector of real numbers for which to find an integer relation.
    # The first term is arctan(1) which corresponds to pi/4.
    terms = [mpmath.atan(1)] + [mpmath.atan(mpmath.mpf(1)/x) for x in denominators]

    # Use the PSLQ algorithm to find a list of integer coefficients 'coeffs'
    # such that coeffs[0]*terms[0] + coeffs[1]*terms[1] + ... = 0
    print("Finding integer relation using the PSLQ algorithm...")
    coeffs = mpmath.pslq(terms)
    print("PSLQ algorithm finished.")
    
    # The relation found is:
    # coeffs[0]*arctan(1) + coeffs[1]*arctan(1/122) + ... = 0
    # The problem's equation is:
    # n*arctan(1) - c1*arctan(1/122) - ... = 0
    #
    # By comparing these, we get:
    # n = coeffs[0]
    # c_i = -coeffs[i+1] for i=1 to 6
    
    n = coeffs[0]
    c = [-k for k in coeffs[1:]]

    # The problem asks for the smallest positive value of n.
    # If PSLQ returned a relation with a negative n, we just flip the signs of all coefficients.
    if n < 0:
        n = -n
        c = [-val for val in c]
    
    c1, c2, c3, c4, c5, c6 = c

    print("\nFound the unique integer coefficients for the smallest positive n:")
    print(f"n = {n}")
    print(f"c1 = {c1}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"c4 = {c4}")
    print(f"c5 = {c5}")
    print(f"c6 = {c6}")

    # Output the full equation with the found numbers.
    print("\nThe resulting equation is:")
    print(f"{n} * pi/4 = {c1}*arctan(1/{denominators[0]}) + {c2}*arctan(1/{denominators[1]}) + {c3}*arctan(1/{denominators[2]}) + {c4}*arctan(1/{denominators[3]}) + {c5}*arctan(1/{denominators[4]}) + {c6}*arctan(1/{denominators[5]})")

    # Final answer in the required format for parsing.
    # Do not print the <<<...>>> part, it is a marker for the final answer.
    # print(f"<<<{n},{c1},{c2},{c3},{c4},{c5},{c6}>>>")

if __name__ == '__main__':
    solve_machin_like_formula()
