import math

def calculate_geodesic_length_bound():
    """
    Calculates the smallest upper bound for the length of a closed geodesic
    on a two-sphere with a given surface area, based on the optimal
    systolic inequality for the 2-sphere.
    """
    # The given surface area of the Riemannian two-sphere.
    surface_area = 8.0

    # According to the optimal systolic inequality by Ivanov, Panov, and Petrunin (2020),
    # for a 2-sphere with area A, there exists a simple closed geodesic of length L such that:
    # L^2 <= (4/pi) * A
    # This gives an upper bound on L: L <= 2 * sqrt(A / pi)

    # Calculate the constant part of the formula
    two = 2
    
    # Calculate the upper bound for the length L.
    length_bound = two * math.sqrt(surface_area / math.pi)

    # Print the final equation with each number and the result, as requested.
    print(f"The upper bound L is calculated by the formula: 2 * sqrt(Area / pi)")
    print(f"Substituting the values: {two} * sqrt({surface_area} / {math.pi}) = {length_bound}")

calculate_geodesic_length_bound()

<<<3.19153824321146>>>