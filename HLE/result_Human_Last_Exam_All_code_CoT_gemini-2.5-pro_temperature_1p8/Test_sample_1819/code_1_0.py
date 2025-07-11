import sys

# On some systems, fractions are not displayed nicely. Use this for better output.
if sys.stdout.encoding != 'UTF-8':
    def pretty_print(num, den):
        return f"{num}/{den}"
else:
    def pretty_print(num, den):
        superscript_map = {
            "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶",
            "7": "⁷", "8": "⁸", "9": "⁹"
        }
        num_str = str(num)
        den_str = str(den)
        return f"{num_str}/{den_str}"

def solve_energy_flow():
    """
    Calculates the energy flow through the yellow sides of the pyramid by
    printing a step-by-step derivation of the surface integral calculation.
    """

    print("To find the energy flow, we must calculate the flux of the vector field F through the two yellow surfaces.")
    print("Flux is given by the surface integral: ∫∫ F ⋅ dS.")

    print("\n--- 1. Define the Geometry and Surfaces ---")
    print("Pyramid Base Vertices: (±1, ±1, 0)")
    print("Pyramid Apex: (0, 0, 4)")
    print("Yellow Sides: Front face (in +x direction) and Back face (in -x direction).")
    print("Vector Field F = (3x³y²z, 3x²y³, z)")

    print("\n--- 2. Plan the Calculation ---")
    print("Due to symmetry, the flux through the front and back faces is identical.")
    print("We will calculate the flux for the front face (Flux₁) and the total flow will be 2 * Flux₁.")

    print("\n--- 3. Flux Calculation for the Front Face ---")
    print("a) The front face is a triangle with vertices (1,-1,0), (1,1,0), and (0,0,4).")
    print("   Its equation is: 4x + z = 4.")

    print("\nb) We parameterize the surface by projecting it onto the yz-plane:")
    print("   x = 1 - z/4")
    print("   The domain for the integration is a triangle in the yz-plane defined by:")
    print("   0 ≤ z ≤ 4  and  -(1 - z/4) ≤ y ≤ (1 - z/4)")

    print("\nc) The outward normal vector element dS is (1, 0, 1/4) dy dz.")

    print("\nd) The integrand F ⋅ dS is calculated on the surface:")
    integrand_str = "3(1 - z/4)³y²z + z/4"
    print(f"   F ⋅ dS = {integrand_str}")

    print("\ne) The double integral for the flux (Flux₁) is:")
    print(f"   Flux₁ = ∫(from z=0 to 4) [ ∫(from y=-(1-z/4) to 1-z/4) ({integrand_str}) dy ] dz")

    print("\nf) First, integrate with respect to y:")
    inner_integral_result = "2z(1 - z/4)⁶ + (z/2)(1 - z/4)"
    print(f"   The result of the inner integral is: {inner_integral_result}")

    print("\ng) Now, integrate the result with respect to z from 0 to 4:")
    print(f"   Flux₁ = ∫(from z=0 to 4) [ {inner_integral_result} ] dz")
    print("   This integral is solved by splitting it into two parts:")
    
    # Values from manual calculation
    part_A_val_num = 4
    part_A_val_den = 7
    part_B_val_num = 4
    part_B_val_den = 3
    
    print(f"   Part A: ∫ 2z(1-z/4)⁶ dz = {part_A_val_num}/{part_A_val_den}")
    print(f"   Part B: ∫ (z/2)(1-z/4) dz = {part_B_val_num}/{part_B_val_den}")

    flux1_num = part_A_val_num * part_B_val_den + part_B_val_num * part_A_val_den
    flux1_den = part_A_val_den * part_B_val_den
    print(f"   Adding the parts: Flux₁ = {part_A_val_num}/{part_A_val_den} + {part_B_val_num}/{part_B_val_den} = {flux1_num}/{flux1_den}")
    
    print("\n--- 4. Final Answer ---")
    total_flux_num = 2 * flux1_num
    total_flux_den = flux1_den
    print("The total energy flow is twice the flux through one yellow face.")
    
    print("\nFinal Equation:")
    print(f"Total Flow = 2 * (Flux₁)")
    print(f"           = 2 * ( {part_A_val_num}/{part_A_val_den} + {part_B_val_num}/{part_B_val_den} )")
    print(f"           = 2 * ( {flux1_num}/{flux1_den} )")
    print(f"           = {total_flux_num}/{total_flux_den}")

solve_energy_flow()
<<<80/21>>>