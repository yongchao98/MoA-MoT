import math

def solve_parrot_math():
    """
    This function formulates the calculation for the parrot and determines the answer code.
    """
    # Step 1: Define the components of the calculation using fractions with small integers.
    # Volume constant from the formula V = (4/3) * pi * r^3
    vol_const_num = 4
    vol_const_den = 3

    # Radius r = 0.5 cm = 1/2 cm
    radius_num = 1
    radius_den = 2

    # Density rho = 0.9 kg/cm^3 = 9/10 kg/cm^3
    density_num = 9
    density_den = 10

    # Step 2: Find a suitable approximation for pi.
    # The error on the final mass must be <= 10%, which means the error on pi must also be <= 10%.
    # 0.9 * pi <= pi_approx <= 1.1 * pi
    # 0.9 * 3.14159... <= pi_approx <= 1.1 * 3.14159...
    # 2.827... <= pi_approx <= 3.456...
    # The whole number 3 (as 3/1) is a simple choice that fits this range.
    pi_approx_num = 3
    pi_approx_den = 1

    # Step 3: Present the final calculation for the parrot.
    # The problem asks to output each number in the final equation.
    # The equation is: mass = (4/3) * pi_approx * r^3 * density
    # mass = (4/3) * 3 * (1/2)^3 * (9/10)
    print("Here is the calculation for the parrot:")
    print(f"mass = ({vol_const_num} / {vol_const_den}) * {pi_approx_num} * ({radius_num} / {radius_den}) * ({radius_num} / {radius_den}) * ({radius_num} / {radius_den}) * ({density_num} / {density_den})")

    # Step 4: Determine the largest integer 'z' used in the calculation.
    all_integers = [
        vol_const_num, vol_const_den,
        pi_approx_num, pi_approx_den,
        radius_num, radius_den,
        density_num, density_den
    ]
    z = max(all_integers)
    # The question format is Yz if possible, N0 otherwise.
    # Since z <= 10, it is possible.
    # The final answer is constructed but not printed here, as per instructions.
    final_answer = f"Y{z}"

solve_parrot_math()