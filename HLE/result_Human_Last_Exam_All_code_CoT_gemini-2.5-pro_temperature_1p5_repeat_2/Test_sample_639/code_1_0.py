import math

def calculate_and_print_distances():
    """
    Calculates and prints the possible normalized distances (r) between spheres
    in a 2D jammed packing for r <= 3.

    The distances are derived from the geometric arrangements of particles in
    stable packings, such as square and hexagonal lattices.
    """
    print("For a 2D system of packed hard spheres, the possible distances between")
    print("any two particles are determined by the geometry of their local packing.")
    print("Considering arrangements from square and hexagonal lattices for r <= 3, the set of distances is:\n")

    # List of tuples, where each tuple contains:
    # (value, list_of_numbers_in_equation, equation_template_string)
    distance_data = [
        (1.0, [1], "{}"),
        (math.sqrt(2), [1, 1], "sqrt({}^2 + {}^2)"),
        (math.sqrt(3), [3], "sqrt({})"),
        (2.0, [1, 1], "{} + {}"),
        (math.sqrt(5), [2, 1], "sqrt({}^2 + {}^2)"),
        (math.sqrt(7), [7], "sqrt({})"),
        (math.sqrt(8), [2, 2], "sqrt({}^2 + {}^2)"),
        (3.0, [1, 1, 1], "{} + {} + {}"),
    ]

    # Iterate through the data and print in the specified format
    for value, nums, template in distance_data:
        # Construct the equation string with the specific numbers
        equation_str = template.format(*nums)
        # Print the final result line
        print(f"r = {equation_str} = {value:.2f}")

calculate_and_print_distances()