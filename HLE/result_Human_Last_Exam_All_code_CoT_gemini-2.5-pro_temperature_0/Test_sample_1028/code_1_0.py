import numpy as np

def demonstrate_overparameterization():
    """
    This function demonstrates why using overly complex models, like high-degree
    polynomials as suggested in option C, does not help with identifiability
    and can lead to poor inference. This is analogous to the problem in the
    birth-death model.
    """
    print("--- Demonstrating Over-parameterization (The problem with option C) ---")

    # 1. Define a simple, true underlying process (a line).
    # This is like the "true" but unknown speciation/extinction rate trend.
    true_slope = 2.5
    true_intercept = 10.0

    # 2. Generate a few data points with some noise.
    # In phylogenetics, the "data points" are the branching times in the tree.
    # We have a limited number of these, just like we have a limited number of points here.
    np.random.seed(0)
    num_points = 6 # A small number of data points
    x_data = np.linspace(0, 5, num_points)
    y_noise = np.random.normal(0, 1.5, num_points)
    y_data = true_slope * x_data + true_intercept + y_noise

    print(f"\nWe have {num_points} noisy data points generated from a simple line.")

    # 3. Fit an overly complex model, as suggested in option C.
    # A polynomial of degree 5 has 6 coefficients. With only 6 data points,
    # the polynomial will fit the noise perfectly, rather than finding the true line.
    degree = 5
    print(f"Fitting a polynomial of degree {degree}, which has {degree + 1} coefficients.")
    
    # This is the core of the fitting process
    coefficients = np.polyfit(x_data, y_data, degree)

    # 4. Print the resulting equation.
    # The coefficients will be highly specific to the random noise in our data.
    # A slightly different dataset would yield a completely different polynomial.
    # This means the parameters are "unidentifiable" - they are not stable.
    print("\nThe model is over-parameterized, fitting the noise perfectly.")
    print("The resulting polynomial equation is: y = ...")
    equation_parts = []
    for i, coeff in enumerate(coefficients):
        power = degree - i
        # We print each number in the equation as requested
        print(f"Coefficient for x^{power}: {coeff:.4f}")
        equation_parts.append(f"({coeff:.4f} * x^{power})")
    
    final_equation = " + ".join(equation_parts).replace("+ (-", "- (")
    print(f"\nFull Equation: y = {final_equation}")
    print("\nThis complex equation fits the 6 data points perfectly but does not represent the true underlying linear trend.")
    print("Similarly, fitting a highly flexible birth-death model (Option C) without additional data (like fossils) or constraints (priors) does not solve the identifiability problem; it makes it worse.")

demonstrate_overparameterization()

<<<C>>>