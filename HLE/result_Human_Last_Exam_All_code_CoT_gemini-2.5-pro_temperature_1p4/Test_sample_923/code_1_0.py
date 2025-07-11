import sys

# The error in the GDscript is on the line: if "texture" in list(json_data):
# GDScript does not have a built-in `list()` function to convert a Dictionary to an Array,
# unlike Python. This is a common point of confusion for developers switching languages.
# The correct way to check for a key in a Dictionary in GDScript is simply: if "texture" in json_data:
# Because `list` is not a defined function, the Godot parser will report that the identifier
# cannot be found in the current scope as a function.
error_message = "Parser Error: The identifier \"list\" isn't a function in the current scope."

print(error_message, file=sys.stdout)