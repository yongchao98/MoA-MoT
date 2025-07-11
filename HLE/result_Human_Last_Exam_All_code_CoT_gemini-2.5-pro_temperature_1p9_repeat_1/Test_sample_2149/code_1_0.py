# The values for Y1, Y2, Y3, and Y4 are deciphered from the clues.
Y1 = 1859
Y4 = 641

# The "Hall" value refers to the year of the Hall-Héroult process invention.
hall_year = 1886

# Clue 3 specifies two reactants. Their simplest indices are 1 and 1.
index1 = 1
index2 = 1

# The phrase "Y1-Hall" implies a subtraction.
y1_hall_term = Y1 - hall_year

# The indices are summed.
indices_term = index1 + index2

# The phrase "Y4 to the Y1-Hall..." is interpreted as an addition of all terms.
result = Y4 + y1_hall_term + indices_term

# Print the final equation with all its components as requested.
print(f"The interpreted equation is: {Y4} + ({Y1} - {hall_year}) + ({index1} + {index2})")
print(f"The final calculated value is: {result}")