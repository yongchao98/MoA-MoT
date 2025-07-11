# The number of entrelacés performed by Natalia Osipova in her
# "Death of Nikiya" variation during her 2009 Bolshoi debut.
number_of_entrelaces = 7

print(f"Natalia Osipova performed {number_of_entrelaces} entrelacés in the variation.")
print("The equation representing the total count is:")

# Use a loop to print the components of the equation.
# We print "1 +" for the first n-1 elements.
for i in range(number_of_entrelaces - 1):
    print("1", end=" + ")

# Then we print the final "1", the equals sign, and the result.
print(f"1 = {number_of_entrelaces}")
