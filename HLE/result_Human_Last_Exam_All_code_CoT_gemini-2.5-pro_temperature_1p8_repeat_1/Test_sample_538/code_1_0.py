import numpy as np

def calculate_reflection_splitting():
    """
    Calculates and prints the number of Bragg reflections for specific plane
    families in a rhombohedrally distorted perovskite (R3m space group)
    using pseudocubic indexing.
    """
    # The rhombohedral distortion occurs along a <111> cubic direction,
    # which becomes the unique axis. We use [1, 1, 1] as the reference.
    unique_axis = np.array([1, 1, 1])
    
    # We analyze families of planes by checking their representative normals.
    # We choose normals that will represent different orientations after distortion.
    families = {
        "{200}": [[2, 0, 0]],
        "{220}": [[2, 2, 0], [2, -2, 0]],
        "{222}": [[2, 2, 2], [2, 2, -2]]
    }

    print("Calculating the number of split Bragg reflections for a rhombohedral (R3m) system:")
    print("-" * 75)
    
    results = {}
    for family_name, representative_normals in families.items():
        # Store the unique dot product values (related to the angle)
        dot_products = set()

        for normal_vec_list in representative_normals:
            normal_vec = np.array(normal_vec_list)
            # We compare the absolute value of the dot product.
            # Different values mean different angles relative to the unique axis.
            dot_product_val = np.abs(np.dot(normal_vec, unique_axis))
            # We round to avoid floating-point inaccuracies
            dot_products.add(round(dot_product_val, 5))

        # The number of unique dot product values equals the number of split peaks
        results[family_name] = len(dot_products)

    # Output the results in the required equation-like format
    print("Final Results:")
    print(f"Number of reflections for the {{200}} family = {results['{200}']}")
    print(f"Number of reflections for the {{220}} family = {results['{220}']}")
    print(f"Number of reflections for the {{222}} family = {results['{222}']}")

calculate_reflection_splitting()