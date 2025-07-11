import numpy as np

def calculate_age(eigenvalues):
    """
    Calculates the age of a group element from its eigenvalues.
    The age is the sum of the fractional parts of the arguments of the eigenvalues, normalized by 2*pi.
    The fractional part is chosen to be in [0, 1).
    """
    age = 0
    for val in eigenvalues:
        # angle() returns value in [-pi, pi]. We scale to [-0.5, 0.5]
        angle_normalized = np.angle(val) / (2 * np.pi)
        # We need the fractional part in [0, 1)
        fractional_part = angle_normalized if angle_normalized >= 0 else angle_normalized + 1
        # Round to handle potential floating point inaccuracies for simple fractions
        age += round(fractional_part, 8)
    return round(age)

def solve_rank_h2c():
    """
    Solves the problem by counting the number of junior (age 1) conjugacy classes
    of the icosahedral group A_5 in its standard 3D representation.
    """
    # Conjugacy classes of A5 are identified by element order and character trace.
    # Eigenvalues are for the standard 3D representation in SL(3, C).
    conjugacy_classes = {
        "Identity (order 1)": {
            "eigenvalues": [1, 1, 1]
        },
        "Rotation by pi (order 2)": {
            "eigenvalues": [1, -1, -1]
        },
        "Rotation by 2pi/3 (order 3)": {
            "eigenvalues": [1, np.exp(2j * np.pi / 3), np.exp(-2j * np.pi / 3)]
        },
        "Rotation by 2pi/5 (order 5, class 1)": {
            "eigenvalues": [1, np.exp(2j * np.pi / 5), np.exp(-2j * np.pi / 5)]
        },
        "Rotation by 4pi/5 (order 5, class 2)": {
            "eigenvalues": [1, np.exp(4j * np.pi / 5), np.exp(-4j * np.pi / 5)]
        }
    }

    print("Calculating the age for each conjugacy class of the icosahedral group A5...")
    
    junior_class_count = 0
    contributions = []
    
    for name, data in conjugacy_classes.items():
        age = calculate_age(data["eigenvalues"])
        print(f"- Class '{name}': Age = {age}")
        if age == 1:
            junior_class_count += 1
            contributions.append("1")
    
    print("\nAccording to the McKay correspondence:")
    print("rank H^2_c(Y, Q) = b4(Y) = b2(Y) = (Number of junior classes)")
    
    final_equation = " + ".join(contributions)
    print(f"\nThe number of junior classes is the sum of contributions from each class with age 1:")
    print(f"Number of junior classes = {final_equation} = {junior_class_count}")
    
    print("\nTherefore, the final rank is:")
    print(f"rank H^2_c(Y, Q) = {junior_class_count}")

solve_rank_h2c()
>>>4