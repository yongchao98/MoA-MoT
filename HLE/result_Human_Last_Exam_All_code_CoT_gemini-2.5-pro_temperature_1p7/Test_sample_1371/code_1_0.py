import math

def solve_seating():
    """
    Calculates the number of ways to arrange the members at the conference lunch.
    """

    # The problem simplifies to arranging two main blocks:
    # 1. SM_Block: 12 scientists and 4 mathematicians, forced to sit together.
    # 2. CE_Block: 5 classicists and 2 ethicists.

    # We calculate arrangements for two primary block orders: Scientists-then-Mathematicians (S-M) and vice-versa (M-S).
    # The calculation is symmetric for both, so we calculate for one and multiply by 2.

    # There are four sub-cases based on the gender of the person at each end of the SM_Block.
    # The number of ways to arrange the SM_Block and CE_Block depends on these genders.
    # F9 = 9!, F5 = 5!
    # Term_FF = (Ways to arrange SM with F ends) * (Ways to arrange CE with F neighbors)
    #         = (48 * F9) * (6 * F5) = 288 * F9 * F5
    # Term_FM = (Ways to arrange SM with F,M ends) * (Ways to arrange CE with F,M neighbors)
    #         = (24 * F9) * (4 * F5) = 96 * F9 * F5
    # Term_MF = (Ways to arrange SM with M,F ends) * (Ways to arrange CE with M,F neighbors)
    #         = (32 * F9) * (4 * F5) = 128 * F9 * F5
    # Term_MM = (Ways to arrange SM with M,M ends) * (Ways to arrange CE with M,M neighbors)
    #         = (16 * F9) * (2 * F5) = 32 * F9 * F5

    # Summing these terms for one order (e.g., S-M)
    # Total_SM_order = (288 + 96 + 128 + 32) * 9! * 5! = 544 * 9! * 5!
    sum_of_coeffs = 288 + 96 + 128 + 32

    # Final Total = 2 * (for S-M and M-S orders) * Total_SM_order
    #             = 2 * 544 * 9! * 5! = 1088 * 9! * 5!

    f9 = math.factorial(9)
    f5 = math.factorial(5)

    c_sm_ff = 48
    c_sm_fm = 24
    c_sm_mf = 32
    c_sm_mm = 16

    c_ce_ff = 6
    c_ce_fm = 4
    c_ce_mf = 4
    c_ce_mm = 2

    final_coeff = 2 * sum_of_coeffs
    total_arrangements = final_coeff * f9 * f5

    print("The final calculation is based on summing arrangements over possible gender pairings at the boundaries of the Scientist/Mathematician block and the Classicist/Ethicist block.")
    print("The total number of arrangements can be expressed by the formula:")
    print(f"Total = 2 * ( ( {c_sm_ff} * 9! ) * ( {c_ce_ff} * 5! ) + ( {c_sm_fm} * 9! ) * ( {c_ce_fm} * 5! ) + ( {c_sm_mf} * 9! ) * ( {c_ce_mf} * 5! ) + ( {c_sm_mm} * 9! ) * ( {c_ce_mm} * 5! ) )")
    print("\nThis simplifies to:")
    print(f"Total = 2 * ( {sum_of_coeffs} * 9! * 5! )")
    print("\nBreaking it down with the factorial values:")
    print(f"Total = {final_coeff} * {f9} * {f5}")
    print("\nThe total number of ways to arrange the table is:")
    print(f"{total_arrangements}")


solve_seating()
<<<474328883200>>>