import sympy

def solve_geometry_ratio():
    """
    This function calculates the ratio of the area of triangle KMN to the area of triangle ABC.

    Given:
    - D, E, F are the points of tangency of the incircle of acute △ABC. △DEF is the contact triangle.
    - r is the inradius of △DEF.
    - R is the inradius of △ABC.
    - △KMN is the orthic triangle of △DEF.

    The derivation steps are:
    1. It can be shown that △KMN is similar to △ABC.
    2. The ratio of their areas is the square of their similarity ratio, which can be expressed using various corresponding lengths (sides, inradii, circumradii).
       S_KMN / S_ABC = (r_in(KMN) / r_in(ABC))^2 = (r_in(KMN) / R)^2.
    3. The inradius of the orthic triangle (KMN) of a parent triangle (DEF) is given by:
       r_in(KMN) = 2 * R_c(DEF) * cos(D')cos(E')cos(F'), where R_c(DEF) is the circumradius of DEF and D', E', F' are its angles.
    4. The circumcircle of the contact triangle DEF is the incircle of ABC. Thus, R_c(DEF) = R.
       So, r_in(KMN) = 2 * R * cos(D')cos(E')cos(F').
    5. Substituting this into the area ratio:
       S_KMN / S_ABC = (2 * R * cos(D')cos(E')cos(F') / R)^2 = (2 * cos(D')cos(E')cos(F'))^2.
    6. A key identity for this configuration is that 2 * cos(D')cos(E')cos(F') = r / R.
    7. Substituting this identity gives the final result:
       S_KMN / S_ABC = (r / R)^2.
    """
    r, R = sympy.symbols('r R')
    
    # The ratio of the areas S_KMN : S_ABC is (r/R)^2
    ratio_expression = (r**2) / (R**2)
    
    # We are asked to output the final equation.
    # To make it clear, we will print the fraction S_KMN / S_ABC = r^2 / R^2
    s_kmn = sympy.Symbol('S_KMN')
    s_abc = sympy.Symbol('S_ABC')

    final_equation = sympy.Eq(s_kmn / s_abc, ratio_expression)

    # We need to print each number in the final equation.
    # The final equation is S_KMN/S_ABC = r^2/R^2. Let's print this representation.
    print(f"The ratio S_KMN : S_ABC is given by the equation:")
    print("S_KMN / S_ABC = r**2 / R**2")
    # For the final answer format, let's represent it symbolically
    # S_KMN / S_ABC = (r^2) / (R^2) has numbers 2 and 2 in the exponents.
    print("Final Equation:")
    print(f"{s_kmn} / {s_abc} = {r}**{2} / {R}**{2}")

solve_geometry_ratio()