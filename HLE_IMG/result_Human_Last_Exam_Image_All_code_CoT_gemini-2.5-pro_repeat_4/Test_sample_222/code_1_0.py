import math

def solve_contour_integral():
    """
    This function calculates the imaginary part of the sum of the integrals as described in the problem.
    """
    
    # Step 1: Define a function to calculate the residue at poles z = -n.
    # The residue of f(z) at z = -n is (z/(z-3/2)) * Res(Gamma, -n)
    # Res(Gamma, -n) = (-1)^n / n!
    # Res(f, -n) = (-n / (-n - 3/2)) * ((-1)^n / n!) = (2n / (2n + 3)) * ((-1)^n / n!)
    def get_residue_at_neg_n(n):
        numerator = 2 * n * ((-1)**n)
        denominator = (2 * n + 3) * math.factorial(n)
        return numerator / denominator

    # Step 2: Calculate the required residues.
    # Residue at z = 3/2 = 1.5
    res_1_5 = (3/4) * math.sqrt(math.pi)
    
    # Residues at z = -1, -2, -3
    res_neg_1 = get_residue_at_neg_n(1)
    res_neg_2 = get_residue_at_neg_n(2)
    res_neg_3 = get_residue_at_neg_n(3)
    
    # Step 3: Determine the sum of integrals using the Residue Theorem.
    # The sum of integrals I = I_C1 + I_C2.
    # I = 2*pi*i * sum[ (Ind_C1(zk) + Ind_C2(zk)) * Res(f, zk) ] over all poles zk.
    #
    # Poles enclosed by C1: 1.5, -1, -2, -3. All with winding number +1.
    # Poles enclosed by C2: -1, -2. Both with winding number +1.
    #
    # Contribution from z=1.5: Ind_C1(1.5)+Ind_C2(1.5) = 1+0 = 1
    # Contribution from z=-1:  Ind_C1(-1)+Ind_C2(-1) = 1+1 = 2
    # Contribution from z=-2:  Ind_C1(-2)+Ind_C2(-2) = 1+1 = 2
    # Contribution from z=-3:  Ind_C1(-3)+Ind_C2(-3) = 1+0 = 1
    #
    # I = 2*pi*i * (1*Res(f, 1.5) + 2*Res(f, -1) + 2*Res(f, -2) + 1*Res(f, -3))
    # The imaginary part of I is 2*pi times the real sum in the parenthesis.
    
    total_residue_sum = res_1_5 + 2 * res_neg_1 + 2 * res_neg_2 + res_neg_3
    
    imaginary_part = 2 * math.pi * total_residue_sum

    # Step 4: Print the calculation step-by-step, including each number.
    print("The imaginary part of the sum of the integrals is Im(I), calculated as:")
    print("Im(I) = 2 * pi * (Res(f, 1.5) + 2*Res(f, -1) + 2*Res(f, -2) + Res(f, -3))")
    print("\nThe individual residue values are:")
    print(f"Res(f, 1.5) = {res_1_5}")
    print(f"Res(f, -1) = {res_neg_1}")
    print(f"Res(f, -2) = {res_neg_2}")
    print(f"Res(f, -3) = {res_neg_3}")

    print("\nSubstituting these values into the equation:")
    print(f"Im(I) = 2 * pi * ({res_1_5} + 2*({res_neg_1}) + 2*({res_neg_2}) + ({res_neg_3}))")
    print(f"Im(I) = 2 * pi * ({res_1_5} + {2*res_neg_1} + {2*res_neg_2} + {res_neg_3})")
    print(f"Im(I) = 2 * pi * ({total_residue_sum})")
    print("\nFinal Answer:")
    print(f"{imaginary_part}")

solve_contour_integral()