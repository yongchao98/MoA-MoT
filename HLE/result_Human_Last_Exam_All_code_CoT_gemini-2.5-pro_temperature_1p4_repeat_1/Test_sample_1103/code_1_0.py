import sys
from cypari2 import pari

def solve_class_number_problem():
    """
    Finds the number of negative fundamental discriminants with class number 48.
    
    This script iterates through negative integers, checks if they are fundamental
    discriminants, and if so, computes their class number. It counts how many
    have a class number of 48.

    This requires the 'cypari2' library, a Python interface for PARI/GP.
    """
    try:
        # Based on known computational results (M. Watkins), the largest absolute
        # value of a discriminant with class number 48 is 630,527. We set the
        # search limit slightly higher to be comprehensive.
        search_limit = 650000

        print(f"Searching for negative fundamental discriminants with class number 48 up to |d| = {search_limit}...")
        
        count_h48 = 0
        discriminants_found = []

        # pari.qfbclassno(d) is highly optimized. We iterate through d = -n.
        for n in range(1, search_limit + 1):
            d = -n
            
            # A fundamental discriminant is an integer d which is the discriminant
            # of a quadratic number field. We use the optimized pari.isfundamental().
            if pari.isfundamental(d):
                # For d < 0, qfbclassno(d) computes the class number of the
                # imaginary quadratic field Q(sqrt(d)).
                class_number = pari.qfbclassno(d)
                
                if class_number == 48:
                    count_h48 += 1
                    discriminants_found.append(d)
        
        print("\n--- Computation Complete ---")
        
        # The final "equation" is the count of the items. We print each item (discriminant)
        # that was counted.
        for d in discriminants_found:
            print(d)

        print("\n--- Final Answer ---")
        # And we print the final resulting number.
        print(f"The total number of negative fundamental discriminants with class number 48 is: {count_h48}")

    except ImportError:
        print("Error: The 'cypari2' library is required to run this script.", file=sys.stderr)
        print("Please install it using: pip install cypari2", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    solve_class_number_problem()
