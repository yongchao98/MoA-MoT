import math

def calculate_sphere_distances():
    """
    This function calculates and prints the set of possible normalized distances (r)
    between two hard spheres on a 2D plane for r <= 3. The distances are derived
    from fundamental geometric arrangements found in jammed packings.
    The instruction "output each number in the final equation" is interpreted as
    showing the components of the calculation (e.g., sqrt(1^2 + 1^2)).
    """

    # We define a list of tuples, where each tuple contains the expression for the distance
    # and a function to evaluate it.
    distance_calculations = [
        ("1.00", lambda: 1.0),
        ("2.00", lambda: 2.0),
        ("3.00", lambda: 3.0),
        ("sqrt( 1^2 + 1^2 )", lambda: math.sqrt(1**2 + 1**2)),
        ("sqrt(3)", lambda: math.sqrt(3)), # From hexagonal packing
        ("sqrt( 2^2 + 1^2 )", lambda: math.sqrt(2**2 + 1**2)),
        ("sqrt(7)", lambda: math.sqrt(7)), # From hexagonal packing
        ("sqrt( 2^2 + 2^2 )", lambda: math.sqrt(2**2 + 2**2)),
    ]

    # Evaluate the distances and store them with their expressions
    results = []
    for expr, func in distance_calculations:
        value = func()
        if value <= 3.0:
            results.append({'value': value, 'expr': expr})

    # Sort the results based on the calculated distance
    sorted_results = sorted(results, key=lambda x: x['value'])

    print("The set of possible normalized distances r <= 3 are:")
    
    # Print the final results in the required format
    for res in sorted_results:
        # For simple integers, just print the value.
        if res['expr'] in ["1.00", "2.00", "3.00"]:
            print(f"{res['expr']}")
        else:
            # For square roots, print the equation and the result.
            print(f"{res['expr']} = {res['value']:.2f}")

calculate_sphere_distances()
<<<[1.00, 1.41, 1.73, 2.00, 2.24, 2.65, 2.83, 3.00]>>>