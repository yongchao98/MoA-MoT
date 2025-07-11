# The following script prints the final equation for the displacement field (D).
# Variable names are used to represent the physical quantities.
# Ctg: Top gate capacitance per area
# Vtg: Top gate voltage
# Cbg: Back gate capacitance per area
# Vbg: Back gate voltage

# Define the symbols for the equation as strings
symbol_D = "D"
symbol_Ctg = "Ctg"
symbol_Vtg = "Vtg"
symbol_Cbg = "Cbg"
symbol_Vbg = "Vbg"

# The instruction "output each number in the final equation" is interpreted as
# outputting each symbolic component of the equation.
print(f"The final equation for the displacement field is:")
print(f"{symbol_D} = ({symbol_Ctg} * {symbol_Vtg}) + ({symbol_Cbg} * {symbol_Vbg})")