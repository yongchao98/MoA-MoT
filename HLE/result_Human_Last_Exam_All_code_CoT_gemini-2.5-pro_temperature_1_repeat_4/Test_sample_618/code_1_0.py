import math

def calculate_x_transformation(t, x, k1, c3, c6, c4, epsilon):
    """
    Calculates the finite Lie group transformation for the x-coordinate.

    This transformation leaves the heat equation with a logarithmic source term,
    u_t = u_{xx} + (k1*ln(u) + k2)u, invariant.

    Args:
        t (float): The initial time coordinate.
        x (float): The initial spatial coordinate.
        k1 (float): The non-zero parameter from the PDE.
        c3 (float): The group parameter associated with time translation.
        c6 (float): The group parameter associated with space translation.
        c4 (float): The group parameter associated with the more complex generator.
        epsilon (float): The continuous parameter of the Lie group.
    """
    if k1 == 0:
        print("Error: This solution is derived for k1 != 0.")
        return

    print("The transformation for the x-coordinate is found by integrating dx'/d(epsilon) = xi(t').")
    print(f"Here, xi(t') = c6 - (2*c4/k1) * exp(k1*t'), and t' = t + c3*epsilon.")
    print("-" * 20)
    
    if c3 == 0:
        print("Case 1: The parameter c3 is zero.")
        print("The transformation equation is: x' = x + epsilon * (c6 - (2*c4/k1) * exp(k1*t))")
        
        xi = c6 - (2 * c4 / k1) * math.exp(k1 * t)
        x_prime = x + epsilon * xi
        
        print("\nSubstituting the provided values:")
        print(f"Equation: x' = {x} + {epsilon} * ({c6} - (2 * {c4} / {k1}) * exp({k1} * {t}))")
        print(f"Intermediate value xi = {xi}")
        print(f"Final Result: x' = {x_prime}")

    else: # c3 != 0
        print("Case 2: The parameter c3 is non-zero.")
        print("The transformation equation is: x' = x + c6*epsilon - (2*c4 / (k1^2 * c3)) * exp(k1*t) * (exp(k1*c3*epsilon) - 1)")
        
        term1 = c6 * epsilon
        exp_factor = math.exp(k1 * t)
        exp_term = math.exp(k1 * c3 * epsilon) - 1
        factor = (2 * c4) / (k1**2 * c3)
        
        x_prime = x + term1 - factor * exp_factor * exp_term
        
        print("\nSubstituting the provided values:")
        print(f"Equation: x' = {x} + {c6}*{epsilon} - (2*{c4} / ({k1}^2 * {c3})) * exp({k1}*{t}) * (exp({k1}*{c3}*{epsilon}) - 1)")
        print(f"Term 1 (c6*epsilon) = {term1}")
        print(f"Term 2 (factor*exp_factor*exp_term) = {factor * exp_factor * exp_term}")
        print(f"Final Result: x' = {x_prime}")

# --- Example Usage ---
print("Running an example for Case 1 (c3=0):")
calculate_x_transformation(t=1.0, x=2.0, k1=0.5, c3=0, c6=3.0, c4=4.0, epsilon=0.1)
print("\n" + "="*40 + "\n")
print("Running an example for Case 2 (c3!=0):")
calculate_x_transformation(t=1.0, x=2.0, k1=0.5, c3=2.0, c6=3.0, c4=4.0, epsilon=0.1)