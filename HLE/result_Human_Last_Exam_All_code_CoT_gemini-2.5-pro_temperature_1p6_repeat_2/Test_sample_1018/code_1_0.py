import math

def compute_simplicial_volume_of_product(g1, g2):
    """
    Computes the simplicial volume of the product of two surfaces Σ_g1 x Σ_g2.

    The computation relies on properties of the fundamental groups of the surfaces.
    """
    
    print(f"Task: Compute the simplicial volume of the product manifold Σ_{g1} x Σ_{g2}.")
    print("-" * 70)
    print("Step 1: Analyze the fundamental groups of the factor surfaces.")
    print(f"The first surface is Σ_{g1} with genus g1 = {g1}.")
    print(f"The second surface is Σ_{g2} with genus g2 = {g2}.")
    
    # A key property of the fundamental group of a surface Σ_g is its amenability.
    # - If g >= 2, the group π₁(Σ_g) is non-amenable.
    # - If g < 2 (i.e., g=0 or g=1), the group π₁(Σ_g) is amenable.

    # We apply a theorem by Gromov regarding the simplicial volume of product manifolds.
    # The theorem states: ||M₁ × M₂|| = 0 if the fundamental group of M₁ or M₂ is non-amenable.

    if g1 >= 2 or g2 >= 2:
        print("\nStep 2: Check the condition for Gromov's Vanishing Theorem.")
        print(f"The condition is that at least one genus must be 2 or greater.")
        print(f"In this case, g1 = {g1} >= 2, so the condition is met.")
        print("This means the fundamental group π₁(Σ_{g1}) is non-amenable.")
        
        print("\nStep 3: Apply the theorem and state the result.")
        print("Because one of the factors has a non-amenable fundamental group, the simplicial")
        print("volume of the product manifold is zero.")
        
        result = 0
        print("\nFinal Equation:")
        print(f"||Σ_{g1} x Σ_{g2}|| = {result}")

    else: # This block handles cases where both g1 < 2 and g2 < 2
        print("\nStep 2: Check the condition for Gromov's Vanishing Theorem.")
        print(f"The condition that a genus is >= 2 is not met (g1={g1}, g2={g2}).")
        print("In this case, both fundamental groups π₁(Σ_{g1}) and π₁(Σ_{g2}) are amenable.")
        
        print("\nStep 3: Apply the theorem for amenable groups and state the result.")
        print("The fundamental group of the product, π₁(Σ_{g1} x Σ_{g2}), is therefore also amenable.")
        print("A different theorem states that any manifold with an amenable fundamental group")
        print("has a simplicial volume of zero.")

        result = 0
        print("\nFinal Equation:")
        print(f"||Σ_{g1} x Σ_{g2}|| = {result}")

# Set the genera from the problem
genus1 = 31
genus2 = 17

compute_simplicial_volume_of_product(genus1, genus2)
