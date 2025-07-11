# This script requires the 'cypari2' library to be installed.
# You can install it by running the following command in your shell:
# pip install cypari2

def find_discriminants_by_class_number():
    """
    Finds and counts negative fundamental discriminants for a given class number.
    """
    try:
        from cypari2 import pari
    except ImportError:
        print("Error: The 'cypari2' library is required but not found.")
        print("Please install it by running: pip install cypari2")
        return

    class_number_to_find = 48
    discriminants_found = []

    # The largest known discriminant for h=48 has an absolute value of 135547.
    # We will search down to -140000 to ensure we find all of them.
    search_limit = -140000

    # Iterate through negative integers within the search limit.
    for d in range(-3, search_limit, -1):
        # Check if d is a fundamental discriminant.
        if pari.isfundamental(d):
            # If it is, compute its class number.
            class_number = pari.qfbclassno(d)
            if class_number == class_number_to_find:
                discriminants_found.append(d)

    count = len(discriminants_found)
    
    print(f"The target class number is: {class_number_to_find}")
    
    # Per the instructions, we output the numbers that contribute to the final answer.
    # In this case, these are the individual discriminants we found.
    print(f"\nThe {count} negative fundamental discriminants with class number {class_number_to_find} are:")
    # We sort them for readability.
    print(sorted(discriminants_found))
    
    # The final "equation" is the statement of the total count.
    print("\nThe final count is:")
    print(f"Total Number = {count}")


if __name__ == "__main__":
    find_discriminants_by_class_number()