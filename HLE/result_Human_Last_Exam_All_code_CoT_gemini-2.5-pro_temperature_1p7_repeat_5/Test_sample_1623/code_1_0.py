# Define the parameters of the identified torus knot T(p,q)
p = 5
q = 4

# The formula for the maximal Thurston-Bennequin number of the T(p,q) torus knot
# is p*q - p - q.
result = p * q - p - q

# Print the explanation and the calculation, showing each number in the equation.
print(f"The knot is identified as the torus knot T({p},{q}).")
print("The maximal Thurston-Bennequin number is calculated using the formula: p * q - p - q.")
print(f"Plugging in the values for p and q:")
print(f"{p} * {q} - {p} - {q} = {result}")
