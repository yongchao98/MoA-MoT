def solve_sum():
    """
    Calculates the sum of 1/n^2 over the specified set S by computing ζ(6) * ζ(8).
    """

    # Known values for the Riemann Zeta function at even integers.
    # ζ(6) = (π^6) / 945
    zeta_6_num = 1
    zeta_6_den = 945
    zeta_6_pi_power = 6

    # ζ(8) = (π^8) / 9450
    zeta_8_num = 1
    zeta_8_den = 9450
    zeta_8_pi_power = 8

    print("The sum is equal to the product of ζ(6) and ζ(8).")
    print("-" * 20)
    print(f"ζ(6) = (π^{zeta_6_pi_power}) / {zeta_6_den}")
    print(f"ζ(8) = (π^{zeta_8_pi_power}) / {zeta_8_den}")
    print("-" * 20)

    # Calculate the product.
    final_num = zeta_6_num * zeta_8_num
    final_den = zeta_6_den * zeta_8_den
    final_pi_power = zeta_6_pi_power + zeta_8_pi_power

    # Print the final equation with all numbers.
    print("The final equation is:")
    print(f"Sum = ( (π^{zeta_6_pi_power}) / {zeta_6_den} ) * ( (π^{zeta_8_pi_power}) / {zeta_8_den} )")
    print(f"Sum = ({final_num} * π^{final_pi_power}) / {final_den}")
    print("-" * 20)

    # Print the final answer in the requested format.
    print("The final answer is:")
    print(f"({final_num}/{final_den}) * π^{final_pi_power}")

solve_sum()