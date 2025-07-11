# This script prints the error message that the Godot 4.3 IDE will display.
# The error occurs in the GDscript on the line: if "texture" in list(json_data):
# GDScript does not have a built-in global function named 'list'.
# The parser will therefore be unable to find the identifier "list" in the current scope.

error_message = "The identifier \"list\" isn't declared in the current scope."
print(error_message)
