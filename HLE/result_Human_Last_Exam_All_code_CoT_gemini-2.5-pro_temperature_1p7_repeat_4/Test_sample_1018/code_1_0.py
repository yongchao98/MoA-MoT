def compute_simplicial_volume_of_product_surface(g1, g2):
    """
    Computes the simplicial volume of the product of two oriented closed surfaces.

    The simplicial volume of a product of two closed, connected, orientable surfaces
    Sigma_g1 x Sigma_g2 is 0 if both genera g1 and g2 are non-zero.
    It is also 0 if at least one of the surfaces is a sphere (genus 0), because
    the simplicial volume of a sphere is 0, and the simplicial volume of a product
    with a factor of volume 0 is 0.
    Thus, for any g1, g2 >= 0, the simplicial volume is 0.

    This function implements this established mathematical result.
    
    Args:
        g1 (int): The genus of the first surface.
        g2 (int): The genus of the second surface.

    Returns:
        int: The simplicial volume of Sigma_g1 x Sigma_g2.
    """
    
    print("Step 1: Define the genera of the two surfaces.")
    print(f"Genus of the first surface, g1 = {g1}")
    print(f"Genus of the second surface, g2 = {g2}")
    print("-" * 20)
    
    print("Step 2: Apply the theorem for the simplicial volume of a product of surfaces.")
    print("Theorem: The simplicial volume ||Sigma_g1 x Sigma_g2|| is 0 for all g1, g2 >= 0.")
    print("We check the conditions for our specific case.")

    # The result is 0 for all non-negative genera g1 and g2.
    # The check is included for pedagogical purposes.
    if g1 >= 0 and g2 >= 0:
        result = 0
        print(f"Since g1={g1} >= 0 and g2={g2} >= 0, the theorem applies.")
    else:
        # This case is not physically meaningful for closed surfaces but included for completeness.
        result = "Invalid genera"
        print("Genera must be non-negative integers.")
        
    print("-" * 20)
    print("Step 3: State the final result.")
    
    if isinstance(result, int):
        print(f"The simplicial volume of \u03A3_{g1} x \u03A3_{g2} is calculated as:")
        # The prompt asks to output each number in the final equation.
        # The "equation" here is the direct result from the theorem.
        print(f"||\u03A3_{g1} x \u03A3_{g2}|| = {result}")

# Define the genera from the problem
genus_1 = 31
genus_2 = 17

# Compute and print the result
compute_simplicial_volume_of_product_surface(genus_1, genus_2)