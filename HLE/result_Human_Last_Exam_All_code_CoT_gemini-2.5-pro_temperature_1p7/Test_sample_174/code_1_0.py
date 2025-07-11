import fractions

def calculate_euler_characteristic():
    """
    Calculates the orbifold Euler characteristic of the quotient stack [U/G]
    of smooth plane quartics by the action of PGL(3).
    """

    print("The goal is to compute the orbifold Euler characteristic of the moduli space of smooth plane quartics, which we denote as chi([U/G]).")
    print("This space is equivalent to the moduli space of non-hyperelliptic genus 3 curves, M_3^nh.")
    print("We use the additivity of Euler characteristic on the moduli space of all genus 3 curves, M_3:")
    print("chi(M_3^nh) = chi(M_3) - chi(H_3), where H_3 is the hyperelliptic locus.\n")

    # === Step 1: Calculate chi(M_3) ===
    g = 3
    # The 6th Bernoulli number, B_6, is 1/42.
    B_6 = fractions.Fraction(1, 42)
    chi_M3_denominator = (2 * g) * (2 * g - 2)
    chi_M3 = B_6 / chi_M3_denominator

    print("Step 1: Calculate chi(M_3) using the Harer-Zagier formula.")
    print(f"Formula: chi(M_g) = B_2g / (2g * (2g - 2)) for g = {g}.")
    print(f"We need B_6 = {B_6.numerator}/{B_6.denominator}.")
    print(f"chi(M_3) = B_6 / (6 * 4) = B_6 / 24.")
    print(f"chi(M_3) = ({B_6.numerator}/{B_6.denominator}) / {chi_M3_denominator} = {chi_M3.numerator}/{chi_M3.denominator}\n")
    
    # === Step 2: State the value for chi(H_3) ===
    # This value is taken from the mathematical literature (e.g., Faber, Bergstrom, Tommasi).
    chi_H3 = fractions.Fraction(-1, 336)
    
    print("Step 2: Use the known value for the Euler characteristic of the hyperelliptic locus, chi(H_3).")
    print(f"From the literature, chi(H_3) = {chi_H3.numerator}/{chi_H3.denominator}.\n")

    # === Step 3: Compute the final result ===
    print("Step 3: Compute chi([U/G]) = chi(M_3) - chi(H_3).")
    
    # Numbers for the equation
    n1_num, n1_den = chi_M3.numerator, chi_M3.denominator
    n2_num, n2_den = chi_H3.numerator, chi_H3.denominator
    
    result = chi_M3 - chi_H3
    
    # Common denominator for display purposes
    common_den = 1008
    n1_com_num = n1_num * (common_den // n1_den)
    n2_com_num = -n2_num * (common_den // n2_den)

    print(f"chi([U/G]) = ({n1_num}/{n1_den}) - ({n2_num}/{n2_den})")
    print(f"           = {n1_com_num}/{common_den} + {n2_com_num}/{common_den}")
    
    sum_num = n1_com_num + n2_com_num
    print(f"           = {sum_num}/{common_den}")
    print(f"           = {result.numerator}/{result.denominator}")

if __name__ == '__main__':
    calculate_euler_characteristic()
