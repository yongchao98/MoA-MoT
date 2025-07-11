# First, you need to install the cypari2 library.
# You can do this by running the following command in your terminal or command prompt:
# pip install cypari2

from cypari2 import Pari

def find_discriminants_with_class_number(target_class_number):
    """
    Finds all negative fundamental discriminants with a given class number.

    This function iterates through negative integers up to a search limit,
    checks if they are fundamental discriminants, calculates their class number,
    and collects those that match the target.
    """
    # Initialize the PARI/GP engine
    pari = Pari()

    # From M. Watkins' 2004 paper, "Class numbers of imaginary quadratic fields",
    # the largest absolute value of a discriminant for class number 48 is 478,443.
    # We set a safe upper search limit of 500,000.
    search_limit = 500000

    print(f"Searching for negative fundamental discriminants with class number {target_class_number} up to D = {-search_limit}...")

    found_discriminants = []

    # Iterate through negative integers D. We can start from -3.
    for D in range(-3, -search_limit - 1, -1):
        # A quick check can slightly speed up the process: fundamental
        # discriminants D must satisfy D % 4 == 0 or D % 4 == 1.
        if D % 4 == 0 or D % 4 == 1:
            # Use PARI's built-in function to check if D is a fundamental discriminant
            if pari.isfundamental(D):
                # If it is, calculate the class number
                class_no = pari.qfbclassno(D)
                if class_no == target_class_number:
                    found_discriminants.append(D)

    return sorted(found_discriminants)

if __name__ == "__main__":
    target_h = 48
    discriminants = find_discriminants_with_class_number(target_h)
    count = len(discriminants)

    print(f"\nFound {count} negative fundamental discriminants with class number {target_h}.")
    print("The discriminants are:")
    print(discriminants)

    # To satisfy the request for an equation representing the final count,
    # we format the output as a sum of 1s.
    if count > 0:
        # We print "each number in the final equation" by representing the count as a sum of 1s.
        equation_str = " + ".join(['1'] * count)
        print(f"\nThe final count is derived from the equation: {equation_str} = {count}")
    else:
        print(f"\nNo discriminants found with class number {target_h} within the search range.")
