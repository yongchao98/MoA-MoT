# Define the Chern numbers of the two insulators
c1 = 1
c2 = 1

# The Chern number of a combined system is the sum of the individual Chern numbers.
c_total = c1 + c2

# Print the final equation and the result
print(f"The Chern number of the first insulator is: {c1}")
print(f"The Chern number of the second insulator is: {c2}")
print("The Chern number of the junction (the combined system) is the sum of the individual Chern numbers.")
print(f"Final Equation: {c1} + {c2} = {c_total}")
print(f"The Chern number of the junction is: {c_total}")