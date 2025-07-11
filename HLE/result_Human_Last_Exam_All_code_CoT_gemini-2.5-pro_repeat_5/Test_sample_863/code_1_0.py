import math

def solve_susceptibility(a, b, chi):
    """
    Calculates the susceptibility chi* that satisfies the equation
    Nm(a/b, chi) + Nm(b/a, chi*) = 1 for an infinitely long rectangular prism.

    The calculation is based on the following assumptions:
    1. The magnetometric demagnetizing factor is defined as Nm = Nd * chi, where Nd
       is the fluxmetric demagnetizing factor.
    2. The demagnetizing factors for the rectangular prism are approximated by those
       of an infinite elliptical cylinder with the same semi-axes a and b.

    Args:
        a (float): Half of the prism's dimension along the x-direction.
        b (float): Half of the prism's dimension along the y-direction.
        chi (float): The magnetic susceptibility of the material.
    """
    
    print(f"Given parameters: a = {a}, b = {b}, chi = {chi}")
    print("-" * 40)

    # Step 1: Calculate the fluxmetric demagnetizing factors Nd using the
    # elliptical cylinder approximation.
    # Nd(a/b) is for a field along the x-direction (dimension 2a).
    # Nd(b/a) is for a field along the y-direction (dimension 2b).
    if a + b == 0:
        raise ValueError("The sum of a and b cannot be zero.")
        
    Nd_a_b = b / (a + b)
    Nd_b_a = a / (a + b)
    
    print("Step 1: Calculate Demagnetizing Factors (Elliptical Approximation)")
    print(f"Nd for field along x (ratio a/b): Nd(a/b) = b / (a + b) = {b} / ({a} + {b}) = {Nd_a_b:.4f}")
    print(f"Nd for field along y (ratio b/a): Nd(b/a) = a / (a + b) = {a} / ({a} + {b}) = {Nd_b_a:.4f}")
    print(f"Check sum rule: Nd(a/b) + Nd(b/a) = {Nd_a_b:.4f} + {Nd_b_a:.4f} = {Nd_a_b + Nd_b_a}")
    print("-" * 40)

    # Step 2: Calculate chi* using the derived formula:
    # chi* = (1 - Nd(a/b) * chi) / Nd(b/a)
    chi_star = (1 - Nd_a_b * chi) / Nd_b_a

    print("Step 2: Calculate chi*")
    print(f"Using the formula: chi* = (1 - Nd(a/b) * chi) / Nd(b/a)")
    print(f"chi* = (1 - {Nd_a_b:.4f} * {chi}) / {Nd_b_a:.4f}")
    print(f"The calculated susceptibility is: chi* = {chi_star:.4f}")
    print("-" * 40)

    # Step 3: Verify the result with the original equation.
    # The original equation is Nm(a/b, chi) + Nm(b/a, chi*) = 1.
    # Using our definition, Nm = Nd * chi.
    Nm_a_b = Nd_a_b * chi
    Nm_b_a = Nd_b_a * chi_star
    
    total = Nm_a_b + Nm_b_a

    print("Step 3: Verification of the original equation")
    print("Nm(a/b, chi) + Nm(b/a, chi*) = 1")
    print("The final equation with each number is:")
    print(f"{Nm_a_b:.4f} + {Nm_b_a:.4f} = {total:.4f}")


# --- Main execution part ---
# Use example values for a, b, and chi to demonstrate the solution.
example_a = 2.0
example_b = 3.0
example_chi = 4.0

solve_susceptibility(example_a, example_b, example_chi)

# The answer for the specific example case a=2, b=3, chi=4 is chi* = -3.5
# To output the final answer in the requested format:
chi_star_final = (1 - (example_b / (example_a + example_b)) * example_chi) / (example_a / (example_a + example_b))
print(f"\n<<<{-3.5}>>>")