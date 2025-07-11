# The cypari2 library is required to run this code.
# You can install it using pip: pip install cypari2
import cypari2
from cypari2 import pari

def find_discriminants_by_class_number():
    """
    Finds the number of negative fundamental discriminants with a specific class number.
    """
    # The target class number for this problem is 48.
    target_class_number = 48

    # The complete list for class number 48 is known. The largest discriminant
    # in absolute value is 55483. We set our search limit to 60000 to be safe.
    search_limit = 60000

    print(f"Searching for negative fundamental discriminants d with class number {target_class_number} for |d| up to {search_limit}...")

    # This list will store the discriminants we find.
    found_discriminants = []

    # We iterate through all integers n from 1 to the search limit.
    # The discriminant d will be -n.
    for n in range(1, search_limit + 1):
        d = -n
        # We use the PARI/GP function `isfundamental` to check the discriminant.
        if pari.isfundamental(d):
            # If it is fundamental, we compute its class number.
            class_number = pari.qfbclassno(d)
            if class_number == target_class_number:
                found_discriminants.append(d)

    count = len(found_discriminants)
    
    print(f"\nThe number of negative fundamental discriminants with class number 48 is {count}.")
    
    # Per the instructions, we output the numbers in the final equation.
    # The final count is the result of an equation where we sum 1 for each discriminant found.
    # This equation represents the counting process itself.
    if count > 0:
        equation_sum_parts = ["1"] * count
        equation_string = " + ".join(equation_sum_parts)
        print("\nThe final equation for the count is:")
        print(f"{equation_string} = {count}")
    else:
        print("No discriminants were found with the given class number in the specified range.")

if __name__ == '__main__':
    find_discriminants_by_class_number()