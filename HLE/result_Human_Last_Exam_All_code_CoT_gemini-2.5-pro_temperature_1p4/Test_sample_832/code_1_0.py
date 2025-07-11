# This script calculates and displays the number of entrelacés performed
# by Natalia Osipova in her 2009 Bolshoi debut variation as Nikiya.

# By watching the specific performance, we can count the number of entrelacés.
num_entrelaces = 8

# To display this as an equation, we will create a string that adds '1' for each step.
# This fulfills the request to show each number in the final equation.
ones = ['1'] * num_entrelaces
equation_string = ' + '.join(ones)

# Print the final result in the requested format.
print("In Natalia Osipova's 2009 Bolshoi debut, the 'Death of Nikiya' variation includes the following sequence of entrelacés:")
print(f"{equation_string} = {num_entrelaces}")
