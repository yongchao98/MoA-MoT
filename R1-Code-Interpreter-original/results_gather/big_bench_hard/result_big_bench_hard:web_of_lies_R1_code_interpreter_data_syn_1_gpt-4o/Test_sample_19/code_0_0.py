# Define the truth values
Fidel = False  # Fidel lies
Sal = not Fidel  # Sal says Fidel tells the truth
Osvaldo = not Sal  # Osvaldo says Sal lies
Alexis = not Osvaldo  # Alexis says Osvaldo lies
Millie = Alexis  # Millie says Alexis tells the truth

# Check if Millie tells the truth
print(Millie)