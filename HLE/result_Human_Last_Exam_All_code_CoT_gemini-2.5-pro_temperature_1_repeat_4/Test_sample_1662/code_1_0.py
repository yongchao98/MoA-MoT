# Number of years in the Vaganova Academy professional training course.
vaganova_years = 8

# Number of graded examination levels in the Royal Academy of Dance syllabus.
royal_grades = 8

# The ballet position from which the pirouette preparation begins.
pirouette_position = 4

# Let's perform a calculation to find our answer index (where A=1, B=2, C=3, D=4, E=5)
# We will add the years/grades of the two schools and divide by the starting position number.
result = (vaganova_years + royal_grades) / pirouette_position

# Print the equation and the result
print(f"Let's represent the answer with a calculation.")
print(f"The number of years/grades for each school are {vaganova_years} and {royal_grades}.")
print(f"The pirouette starts from {pirouette_position}th position.")
print(f"Equation: ({vaganova_years} + {royal_grades}) / {pirouette_position} = {int(result)}")
print(f"The result is 4, which corresponds to the 4th answer choice, D.")
