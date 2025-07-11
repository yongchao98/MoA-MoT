# The artifact is a Ptolemaic intaglio, with features characteristic of
# the reign of Ptolemy II Philadelphus (285–246 B.C.). We can use a
# representative year from this period to calculate the century.
representative_year_bc = 250

# The equation to determine the century from a B.C. year is:
# (year // 100) + 1
# Here are the numbers from that equation:
number_1 = representative_year_bc
number_2 = 100
number_3 = 1

# Perform the calculation
century = (number_1 // number_2) + number_3

# The calculated century is 3. Now we will format and print the answer.
# The following line prints the number from our equation as part of the final answer.
print(f"{century}rd century B.C.")