import cypari2
from cypari2 import pari

def count_discriminants_by_class_number(target_class_number):
    """
    Counts the number of negative fundamental discriminants with a given class number.

    This function relies on the 'cypari2' library, a Python wrapper for the
    PARI/GP number theory library.

    Args:
        target_class_number (int): The class number to search for.

    Returns:
        int: The number of discriminants found.
    """
    print(f"Searching for negative fundamental discriminants with class number {target_class_number}...")
    
    # According to M. Watkins' 2004 paper, the largest |d| for h=48 is 243,831.
    # We set our search limit based on this known result.
    if target_class_number == 48:
        limit = 243832
    else:
        # For other class numbers, a limit would need to be determined.
        # This is a non-trivial problem. We'll set a default for demonstration.
        # This code is specifically tailored for class number 48.
        print("Warning: Search limit is optimized for class number 48.")
        limit = 300000 

    count = 0
    found_discriminants = []

    # Iterate through the absolute values of the discriminants
    for D in range(3, limit):
        d = -D
        
        # Optimization: Discriminants must be congruent to 0 or 1 mod 4.
        if d % 4 != 0 and d % 4 != 1:
            continue

        # Use PARI/GP to check if d is a fundamental discriminant
        if pari.isfundamental(d):
            # If so, calculate its class number
            class_number = pari.qfbclassno(d)
            
            if class_number == target_class_number:
                count += 1
                found_discriminants.append(d)

    # The user did not want the list of discriminants, but it's good practice
    # to show them for verification. We can keep it commented out for the final run.
    # print(f"Found discriminants: {found_discriminants}")
                
    return count

# Set the target class number
TARGET_H = 48

# Run the calculation
final_count = count_discriminants_by_class_number(TARGET_H)

# The final answer equation as requested
print(f"The number of negative fundamental discriminants with class number 48 is: {final_count}")
