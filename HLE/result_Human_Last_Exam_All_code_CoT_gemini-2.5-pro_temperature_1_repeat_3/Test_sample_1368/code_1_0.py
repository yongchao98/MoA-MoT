# The plan is to construct two perpendicular diameters.
# The endpoints of these diameters will be the vertices of the inscribed square.

# 1. 'L': Draw a Line through the center O and the given point A on the circumference.
#    This line defines the first diameter and creates the opposite point C.
step1 = "L"

# 2. 'C': Draw a Circle centered at point A, with the radius defined by point C.
step2 = "C"

# 3. 'C': Draw another Circle centered at point C, with the radius defined by point A.
#    The intersection of these two circles defines the line for the perpendicular diameter.
step3 = "C"

# 4. 'L': Draw a Line connecting the two intersection points from the previous steps.
#    This line is the second diameter, perpendicular to the first. Its endpoints on the
#    original circle are the other two vertices of the square.
step4 = "L"

# The final sequence of commands is the concatenation of these steps.
final_sequence = step1 + step2 + step3 + step4

# Print the resulting shortest sequence.
print(final_sequence)