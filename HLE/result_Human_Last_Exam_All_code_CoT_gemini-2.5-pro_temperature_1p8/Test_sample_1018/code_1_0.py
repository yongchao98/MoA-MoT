def compute_simplicial_volume_of_product():
    """
    Computes the simplicial volume of the product of two oriented closed surfaces.

    The simplicial volume of a single oriented closed surface of genus g >= 1
    is ||Σ_g|| = 4g - 4.

    However, a key theorem in topology states that the simplicial volume of
    the product of two such surfaces, Σ_g x Σ_h, is 0 if g >= 1 and h >= 1.
    """
    
    # Genus of the first surface
    g = 31
    
    # Genus of the second surface
    h = 17

    # Check if the conditions for the theorem are met
    if g >= 1 and h >= 1:
        # According to the theorem, the simplicial volume of the product is 0.
        result = 0
    else:
        # If one of the surfaces is a sphere (genus 0), its simplicial volume is 0,
        # and the product volume is also 0.
        result = 0

    # Print the final result as an equation including the input numbers
    print(f"Let Σ_g and Σ_h be oriented closed surfaces of genus g and h.")
    print(f"For g = {g} and h = {h}, we compute the simplicial volume of their product, ||Σ_{g} x Σ_{h}||.")
    print(f"According to a theorem for products of surfaces (where g,h >= 1), the result is:")
    print(f"||Σ_{g} x Σ_{h}| = {result}")

compute_simplicial_volume_of_product()