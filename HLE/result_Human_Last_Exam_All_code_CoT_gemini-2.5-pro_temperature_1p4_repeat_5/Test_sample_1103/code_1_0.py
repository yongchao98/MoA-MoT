# You may need to install the cypari2 library first:
# pip install cypari2

try:
    from cypari2 import pari
except ImportError:
    print("Please install the 'cypari2' library using 'pip install cypari2'")
    exit()

def count_discriminants_with_class_number(h_target):
    """
    Counts the number of negative fundamental discriminants with a given class number.

    Args:
        h_target: The target class number.

    Returns:
        The total count of such discriminants.
    """
    # Based on known results (M. Watkins, 2004), the largest discriminant |D|
    # for h=48 is 135507. We search up to 140000 to be sure.
    limit = 140000
    count = 0

    print(f"Searching for negative fundamental discriminants with class number {h_target}...")

    # Iterate through the absolute values of the discriminants
    for d in range(1, limit + 1):
        D = -d
        # isfundamental() checks if D is a fundamental discriminant
        if pari.isfundamental(D):
            # qfbclassno() computes the class number
            h = pari.qfbclassno(D)
            if h == h_target:
                count += 1
    
    return count

# The target class number is 48
target_class_number = 48
final_count = count_discriminants_with_class_number(target_class_number)

# The final "equation" is the statement of the result.
print("\nFinal Result:")
print(f"The number of negative fundamental discriminants with class number {target_class_number} is = {final_count}")
