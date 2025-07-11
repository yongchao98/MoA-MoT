import sympy as sp

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth moment to the average of the dimensionless
    conductance for disordered Majorana wires at the critical point.
    """
    # Step 1: Define the symbolic variable for conductance 'g'
    g = sp.Symbol('g', real=True, positive=True)

    # Step 2: Define the arcsine probability distribution P(g)
    # P(g) = 1 / (pi * sqrt(g * (1 - g)))
    P_g = 1 / (sp.pi * sp.sqrt(g * (1 - g)))

    # Step 3: Calculate the average conductance <g> (the first moment)
    # This is the integral of g * P(g) from 0 to 1.
    avg_g_integrand = g * P_g
    avg_g = sp.integrate(avg_g_integrand, (g, 0, 1))

    # Step 4: Calculate the fourth moment of conductance <g^4>
    # This is the integral of g^4 * P(g) from 0 to 1.
    fourth_moment_g_integrand = g**4 * P_g
    fourth_moment_g = sp.integrate(fourth_moment_g_integrand, (g, 0, 1))

    # Step 5: Calculate the final ratio
    if avg_g == 0:
        ratio = sp.oo # Handle division by zero, though not expected here
    else:
        ratio = fourth_moment_g / avg_g

    # Step 6: Print the results clearly as requested
    print("The dimensionless conductance 'g' for a disordered Majorana wire at the critical point follows the arcsine distribution:")
    print("P(g) = 1 / (pi * sqrt(g * (1 - g)))")
    print("-" * 50)
    print(f"The average value of the conductance, <g>, is calculated as:")
    print(f"<g> = ∫[0,1] g * P(g) dg = {avg_g}")
    print("-" * 50)
    print(f"The fourth statistical moment, <g^4>, is calculated as:")
    print(f"<g^4> = ∫[0,1] g^4 * P(g) dg = {fourth_moment_g}")
    print("-" * 50)
    print("The ratio between the fourth moment and the average value is therefore:")
    # The final equation with all the numbers is printed here
    print(f"<g^4> / <g> = ({fourth_moment_g}) / ({avg_g}) = {ratio}")
    print(f"\nThe numerical value of the ratio is {ratio.evalf()}.")

if __name__ == '__main__':
    solve_conductance_ratio()