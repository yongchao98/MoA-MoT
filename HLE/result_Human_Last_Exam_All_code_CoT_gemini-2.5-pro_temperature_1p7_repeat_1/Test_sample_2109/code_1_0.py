import sympy
from sympy import Symbol, integrate, pi, Rational

def solve_problem():
    """
    Solves the entire problem step-by-step using symbolic math.
    """
    x = Symbol('x')

    # Step 1 & 2: Define the functions T1, T2, and Phi based on the analysis.
    # From the functional equation for T1, we find:
    # T1(y) = (y**2 * (20 - 6*y**2)) / (3 * (y**2 + 2)**2)
    # We need T1(sqrt(2)*x):
    T1_sqrt2x = (2*x**2 * (5 - 3*x**2)) / (3 * (x**2 + 1)**2)

    # From the functional equation for T2, assuming a polynomial solution leads
    # to a linear function T2(x) = -x + 2 (by correcting a likely typo in the problem statement).
    T2 = -x + 2

    # From the fractional derivative constraints, we find Phi(z) = z/2 + 1.
    # Phi_func = lambda z: Rational(1, 2) * z + 1

    # Step 3: Formulate the total energy integral.
    # E_total = integral_0^1 [ Phi(T1(sqrt(2)x) + T2(x)) ] dx
    # E_total = integral_0^1 [ 1/2 * (T1(sqrt(2)x) + T2(x)) + 1 ] dx
    integrand = Rational(1, 2) * (T1_sqrt2x + T2) + 1

    # Step 4: Calculate the integral.
    # It can be split into three parts:
    # 1/2 * integral(T1_sqrt2x) + 1/2 * integral(T2) + integral(1)

    int1 = integrate(T1_sqrt2x, (x, 0, 1))
    int2 = integrate(T2, (x, 0, 1))
    int3 = integrate(1, (x, 0, 1))

    coeff1 = Rational(1, 2)
    coeff2 = Rational(1, 2)
    
    # Calculate the final total energy
    E_total = coeff1 * int1 + coeff2 * int2 + int3
    
    # Print the detailed breakdown of the calculation
    print("Step 1: The temperature distribution functions are determined as:")
    print(f"T1(sqrt(2)*x) = {T1_sqrt2x}")
    print(f"T2(x) = {T2}")
    print("\nStep 2: The heat flux function is determined as:")
    print("Phi(z) = z/2 + 1")
    print("\nStep 3: The total energy is the integral from 0 to 1 of Phi(T1(sqrt(2)*x) + T2(x)).")
    print("E_total = (1/2) * integral(T1) + (1/2) * integral(T2) + integral(1)")
    print("\nStep 4: Calculating the definite integrals:")
    print(f"Integral of T1 component: {int1}")
    print(f"Integral of T2 component: {int2}")
    print(f"Integral of constant component: {int3}")
    print("\nFinal calculation:")
    print(f"E_total = {coeff1} * ({int1}) + {coeff2} * ({int2}) + {int3}")

    final_symbolic_value = E_total
    final_numeric_value = E_total.evalf()

    print(f"\nThe exact minimum value of the total heat energy is: {final_symbolic_value}")
    print(f"The numerical value is approximately: {final_numeric_value}")

    return final_numeric_value

if __name__ == '__main__':
    result = solve_problem()
    # The final numerical result should be presented clearly at the end.
    # The value is (1 + 7*pi)/12
    # <<<1.91554541579308>>>
    # Let's print the result in the requested format
    final_answer = (1 + 7*pi)/12
    print(f"\nFinal Answer:\n<<<{final_answer.evalf()}>>>")