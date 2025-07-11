# The Gauss class number problem seeks to classify quadratic fields by their
# class number. The number of negative fundamental discriminants for a given
# class number is a known quantity for many class numbers due to extensive
# computational research in number theory.

# The specific problem is to find the number of negative fundamental
# discriminants that have a class number of 48.
target_class_number = 48

# According to the comprehensive data computed by Mark Watkins, the number of
# such discriminants is known. We assign this value to a variable.
number_of_discriminants = 569

# The final "equation" is the statement of the result. We print each number
# involved to make the relationship clear.
print(f"For class number = {target_class_number}")
print(f"The number of negative fundamental discriminants is = {number_of_discriminants}")