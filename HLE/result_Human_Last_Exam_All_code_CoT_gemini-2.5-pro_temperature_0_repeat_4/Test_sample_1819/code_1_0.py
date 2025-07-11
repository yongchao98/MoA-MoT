import fractions

def solve_flux_problem():
    """
    Solves the energy flow problem through the yellow sides of a pyramid
    by applying the Divergence Theorem and principles of symmetry.
    """
    print("Step-by-step solution to find the energy flow through the yellow sides of the pyramid.")
    print("-" * 80)

    # Step 1: Define the problem and state the approach using the Divergence Theorem.
    print("1. The Problem:")
    print("   - Vector field F = (3x³y²z, 3x²y³z, z)")
    print("   - A square pyramid with base side 2, height 4, centered at the origin.")
    print("   - Find the flux (energy flow) through two interspersed (yellow) sides.")
    print("\n2. The Approach:")
    print("   - Use the Divergence Theorem: Flux_total = ∫∫∫_V (div F) dV")
    print("-" * 80)

    # Step 2: Calculate the divergence of F.
    print("3. Calculate the Divergence of F:")
    print("   div F = ∂/∂x(3x³y²z) + ∂/∂y(3x²y³z) + ∂/∂z(z)")
    print("         = 9x²y²z + 9x²y²z + 1")
    div_F_str = "18x²y²z + 1"
    print(f"   div F = {div_F_str}")
    print("-" * 80)

    # Step 3: Calculate the volume integral of the divergence to find the total flux.
    print("4. Calculate the Volume Integral of div F (Total Flux):")
    print("   Flux_total = ∫∫∫_V (18x²y²z + 1) dV = ∫∫∫_V 18x²y²z dV + ∫∫∫_V 1 dV")
    
    # Part A: The integral of 1 dV is the volume of the pyramid.
    base_side = 2
    height = 4
    volume = fractions.Fraction(1, 3) * (base_side**2) * height
    print(f"\n   Part A: ∫∫∫_V 1 dV = Volume of the pyramid")
    print(f"   Volume = (1/3) * Base Area * Height = (1/3) * ({base_side}*{base_side}) * {height} = {volume}")

    # Part B: The integral of 18x²y²z dV over the pyramid's volume.
    # This integral has been pre-calculated.
    integral_part_b = fractions.Fraction(16, 7)
    print(f"\n   Part B: The integral of 18x²y²z dV over the pyramid's volume evaluates to 16/7.")
    print(f"   Result = {integral_part_b}")

    # The total flux is the sum of the two parts.
    total_flux = integral_part_b + volume
    print(f"\n   Total Flux = (Result of Part B) + (Result of Part A)")
    print(f"   Flux_total = {integral_part_b} + {volume} = {total_flux}")
    print("-" * 80)

    # Step 4: Calculate the flux through the base of the pyramid.
    print("5. Calculate Flux through the Base:")
    print("   The base lies on the z=0 plane. On this plane, F = (0, 0, 0).")
    flux_base = 0
    print(f"   Therefore, the flux through the base is {flux_base}.")
    print("-" * 80)

    # Step 5: Calculate the total flux through the four slanted sides.
    print("6. Calculate Flux through the 4 Slanted Sides:")
    flux_sides = total_flux - flux_base
    print(f"   Flux_sides = Flux_total - Flux_base = {total_flux} - {flux_base} = {flux_sides}")
    print("-" * 80)

    # Step 6: Calculate the flux through a single slanted side using symmetry.
    print("7. Calculate Flux through a Single Slanted Side:")
    print("   Due to symmetry, the flux is the same for all 4 slanted sides.")
    flux_one_side = flux_sides / 4
    print(f"   Flux_one_side = Flux_sides / 4 = {flux_sides} / 4 = {flux_one_side}")
    print("-" * 80)

    # Step 7: Calculate the final answer for the two yellow sides.
    print("8. Final Calculation: Flux through the Two Yellow Sides:")
    flux_yellow = 2 * flux_one_side
    print(f"   The total energy flow through the two yellow sides is twice the flux through a single side.")
    print(f"   Final Equation: Energy Flow = 2 * ( (Total Flux - Flux through Base) / 4 )")
    print(f"   Energy Flow = 2 * ( ({total_flux} - {flux_base}) / 4 )")
    print(f"               = 2 * ( {flux_sides} / 4 )")
    print(f"               = 2 * {flux_one_side}")
    print(f"               = {flux_yellow}")
    print("-" * 80)
    
    final_answer_val = float(flux_yellow)
    print(f"The final answer is {flux_yellow} (approximately {final_answer_val:.4f}).")

if __name__ == '__main__':
    solve_flux_problem()