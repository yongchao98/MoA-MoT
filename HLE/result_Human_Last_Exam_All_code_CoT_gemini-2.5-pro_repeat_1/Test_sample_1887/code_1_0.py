# This script presents the final answer for the order type of the set X.
# Based on the derivation, the order type is the ordinal ω_{ω+5}.

# The components of this expression are the symbol 'ω' and the number 5.
omega_symbol = "ω"
number_five = 5

# Construct the subscript part of the expression.
subscript = f"{omega_symbol}+{number_five}"

# Construct the full expression for the ordinal.
# The format is base_(subscript).
order_type_expression = f"{omega_symbol}_({subscript})"

# Print the final result.
print("The order type of X is given by the following ordinal expression:")
print(order_type_expression)