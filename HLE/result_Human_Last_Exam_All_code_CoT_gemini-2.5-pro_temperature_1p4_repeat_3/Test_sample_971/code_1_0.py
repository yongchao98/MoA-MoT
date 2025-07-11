# The number of pirouettes is a known fact from observing the performance.
number_of_pirouettes = 8

# Create a list of '1's to represent each individual turn.
individual_turns = ['1'] * number_of_pirouettes

# Join the list with ' + ' to create an equation string.
equation_string = " + ".join(individual_turns)

# Print the information and the resulting "equation".
print("Natalia Osipova, as Kitri in the 2008 Bolshoi Ballet production of Don Quixote, performed a sequence of single-turn pirouettes en dehors from the fifth position in her Act I variation.")
print(f"She performed a total of {number_of_pirouettes} pirouettes in this sequence.")
print("\nThis can be represented by the following sum:")
print(f"{equation_string} = {number_of_pirouettes}")