import sympy

def derive_critical_speed():
    """
    Derives and prints the formula for the critical speed of an oversteering vehicle.
    """
    # 1. Define the symbolic parameters of the single-track model
    a, b, c_f, c_r, m, I = sympy.symbols('a b c_f c_r m I')

    # 2. Define the oversteer condition parameter
    # The vehicle is oversteering if (a * c_f - b * c_r) > 0.
    # This term appears in the denominator of the critical speed formula.
    oversteer_term = a * c_f - b * c_r

    # 3. Construct the numerator of the v_crit^2 formula
    # This comes from the vehicle's geometry and tire properties. L = (a + b) is the wheelbase.
    numerator = c_f * c_r * (a + b)**2

    # 4. Construct the denominator of the v_crit^2 formula
    # This includes the vehicle mass and the oversteer term.
    denominator = m * oversteer_term

    # 5. Form the expression for the square of the critical speed
    v_crit_squared = numerator / denominator

    # 6. The critical speed is the square root of this expression
    v_crit = sympy.sqrt(v_crit_squared)

    # 7. Print the final formula and its components clearly
    print("The formula for the critical speed (v_crit) at which an oversteering vehicle becomes unstable is derived as follows:")
    print("\nFinal Equation:")
    # The user requested each part of the equation be shown. We format the string to look like an equation.
    print(f"v_crit = sqrt( ({numerator}) / ({denominator}) )")

    print("\nWhere the parameters are:")
    print(f"v_crit: Critical speed")
    print(f"a: Distance from Center of Gravity (CG) to the front axle")
    print(f"b: Distance from CG to the rear axle")
    print(f"c_f: Cornering stiffness of the front axle")
    print(f"c_r: Cornering stiffness of the rear axle")
    print(f"m: Vehicle mass")
    # The moment of inertia 'I' cancels out during the derivation of the final formula for v_crit,
    # so it does not appear in the final equation.

if __name__ == '__main__':
    derive_critical_speed()
    # For a numerical example, if we had values:
    # a = 1.2 m, b = 1.6 m, c_f = 90000 N/rad, c_r = 80000 N/rad, m = 1500 kg
    # a*c_f - b*c_r = 1.2*90000 - 1.6*80000 = 108000 - 128000 = -20000.
    # This example vehicle is understeering, so it does not have a critical speed.
    # Let's swap the stiffnesses to make it oversteering:
    # c_f = 80000, c_r = 90000
    # a*c_f - b*c_r = 1.2*80000 - 1.6*90000 = 96000 - 144000 = -48000. Still understeering.
    # Let's try different geometry:
    # a=1.6, b=1.2, c_f = 90000, c_r = 80000
    # a*c_f - b*c_r = 1.6*90000 - 1.2*80000 = 144000 - 96000 = 48000. This is oversteering.
    # num = 90000 * 80000 * (1.6+1.2)**2 = 7.2e9 * 2.8**2 = 5.6448e10
    # den = 1500 * 48000 = 7.2e7
    # v_crit_sq = 5.6448e10 / 7.2e7 = 784
    # v_crit = sqrt(784) = 28 m/s (or about 100.8 km/h)
    
    final_expression = sympy.sqrt((sympy.Symbol('c_f') * sympy.Symbol('c_r') * (sympy.Symbol('a') + sympy.Symbol('b'))**2) / (sympy.Symbol('m') * (sympy.Symbol('a') * sympy.Symbol('c_f') - sympy.Symbol('b') * sympy.Symbol('c_r'))))
    <<<print(final_expression)>>>