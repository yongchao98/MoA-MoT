import sys

def solve_titan_gravity():
    """
    This program simulates the calculation of gravitational force on a probe
    near a black hole, adhering to the constraints of the Titan 6-bit computer.
    It will determine if the calculation is feasible.
    """

    # The maximum value for a 6-bit unsigned integer is 2^6 - 1 = 63.
    MAX_VAL = 63

    print("Titan 6-bit Superconducting Computer Feasibility Test")
    print("------------------------------------------------------")
    print(f"Constraint: All numerators and denominators must be <= {MAX_VAL}.\n")

    # This function simulates a multiplication operation on Titan.
    # It checks if the resulting numerator would exceed the 6-bit limit.
    def titan_multiply(num1, num2, operation_name):
        res_num = num1 * num2
        print(f"Operation: {operation_name}")
        print(f"Attempting to multiply numerators: {num1} * {num2} = {res_num}")
        if res_num > MAX_VAL:
            print(f"-> ERROR: Resulting numerator {res_num} exceeds the 6-bit limit of {MAX_VAL}.")
            return None
        print(f"-> Success: {res_num} is a valid 6-bit integer.")
        return res_num

    print("Objective: Calculate the mass of Pandora, M = (4/3) * pi * rho * R^3.")
    print("We must use the simplest integer approximations to have any chance of success.\n")

    # Approximations for the calculation of Mass (M)
    # Each is represented as a fraction (numerator / denominator).
    # We only focus on the numerators as they are the source of overflow.
    term_4_3_num = 4
    # pi is approx 3.14. Simplest integer approximation is 3.
    pi_num = 3
    # rho is 1.2e3 kg/m^3. Simplest integer approximation is 1e3 kg/m^3.
    rho_num = 1
    # R is 2e6 m, so R^3 is 8e18 m^3.
    R_cubed_num = 8

    print(f"Component numerators for Mass (M):")
    print(f"  - 4/3      -> {term_4_3_num}")
    print(f"  - pi       -> {pi_num} (from approximation pi ≈ 3/1)")
    print(f"  - rho      -> {rho_num} (from approximation rho ≈ 1e3 kg/m^3)")
    print(f"  - R^3      -> {R_cubed_num} (from R = 2e6 m)")
    print("\nThe final numerator for M will be the product of these four numbers.")
    print("-" * 25)

    try:
        # Step 1: Multiply 4 * pi_num
        current_num = titan_multiply(term_4_3_num, pi_num, "4 * pi")
        if current_num is None:
            raise ValueError("Overflow during first multiplication step.")
        print("-" * 25)

        # Step 2: Multiply the result by rho_num
        current_num = titan_multiply(current_num, rho_num, "(4*pi) * rho")
        if current_num is None:
            raise ValueError("Overflow during second multiplication step.")
        print("-" * 25)

        # Step 3: Multiply the result by R_cubed_num
        current_num = titan_multiply(current_num, R_cubed_num, "(4*pi*rho) * R^3")
        if current_num is None:
            raise ValueError("Overflow during final multiplication step for Mass.")
        print("-" * 25)

        # If we somehow succeeded, we would print the result.
        # But the logic shows this is not possible.

    except ValueError:
        print("\n[CONCLUSION]")
        print("The calculation of the planet's mass (M) is not possible.")
        print("Even when using the smallest possible integer approximations for the constants,")
        print("the product of the numerators (4 * 3 * 1 * 8 = 96) exceeds the 6-bit limit of 63.")
        print("Since M cannot be calculated, the subsequent values for the Schwarzschild Radius (Rs)")
        print("and the gravitational force (F) cannot be determined.")
        print("\nThe problem is therefore unsolvable on the Titan architecture.")

    # Final Answer Block
    print("<<<N0>>>")

solve_titan_gravity()