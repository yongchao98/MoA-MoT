import cypari2

# Initialize the PARI/GP environment
pari = cypari2.Pari()

# The class number we are searching for
target_class_number = 48

# Based on published results by M. Watkins, the largest absolute value of a
# fundamental discriminant with class number 48 is 239,299.
# We set our search limit slightly higher to be comprehensive.
search_limit = 240000

# Counter for the number of discriminants found
discriminant_count = 0

# Iterate through all negative integers within our search limit
for d_abs in range(1, search_limit + 1):
    d = -d_abs
    
    # A fundamental discriminant is an integer d which is the discriminant
    # of a quadratic number field. The pari.isfundamental(d) function checks this condition.
    if pari.isfundamental(d):
        
        # If d is a fundamental discriminant, compute its class number.
        # The pari.qfbclassno(d) function calculates the class number for the discriminant d.
        class_number = pari.qfbclassno(d)
        
        # Check if the computed class number matches our target
        if class_number == target_class_number:
            discriminant_count += 1

# Print the final result in the format of an equation,
# showing the input parameter and the computed result.
print(f"Number of negative fundamental discriminants with class number {target_class_number} = {discriminant_count}")