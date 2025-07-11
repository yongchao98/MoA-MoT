import numpy as np

def solve_model():
    """
    This function constructs a parsimonious regression model from the given data,
    identifies significant parameters, and prints the final model equation.
    """
    # 1. Define the data from the problem description
    x1 = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2 = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3 = np.array([-1, -1, -1, -1, 1, 1, 1, 1])
    y = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    n = len(y)

    # 2. Define all possible terms in the full model
    terms = {
        'Intercept': np.ones(n),
        'x1': x1,
        'x2': x2,
        'x3': x3,
        'x1*x2': x1 * x2,
        'x1*x3': x1 * x3,
        'x2*x3': x2 * x3,
        'x1*x2*x3': x1 * x2 * x3
    }
    
    # 3. Calculate the coefficients for the full model using least squares
    # For an orthogonal design, coeff_i = (X_i' * y) / (X_i' * X_i)
    # where X_i is the column for term i, and X_i'*X_i = n = 8 for all terms.
    coefficients = {name: np.sum(term_vector * y) / n for name, term_vector in terms.items()}

    # 4. Select terms for the parsimonious model.
    # We observe the coefficients:
    # b0=51.4, b1=20.8, b2=-0.0, b3=-12.6, b12=0.8, b13=0.1, b23=0.3, b123=10.2
    # There is a clear gap in magnitude. Terms with coefficients |b| > 1.0 seem significant.
    # These are x1, x3, and x1*x2*x3. The intercept is always included.
    
    final_params = {}
    final_params['Intercept'] = coefficients['Intercept']
    for name, coeff in coefficients.items():
        if name != 'Intercept' and abs(coeff) > 1.0:
            final_params[name] = coeff

    # 5. Construct and print the final model equation
    # Start with the intercept
    equation_parts = [f"y = {final_params.pop('Intercept'):.1f}"]

    # Add the other significant terms
    # Sorting for a consistent, readable order (by number of variables in term)
    sorted_terms = sorted(final_params.keys(), key=lambda k: len(k.replace('*','')))
    
    for term in sorted_terms:
        coeff = final_params[term]
        # Format the term with its sign and coefficient
        sign = "-" if coeff < 0 else "+"
        # The term name for x1 is just 'x1', not '1*x1'
        term_str = term
        
        # We need to output each number in the final equation
        # The format string handles the rounding and formatting
        equation_parts.append(f"{sign} {abs(coeff):.1f}*{term_str}")

    # Join all parts with a space for the final equation
    final_equation = " ".join(equation_parts)
    print(final_equation)

solve_model()
<<<y = 51.4 + 20.8*x1 - 12.6*x3 + 10.2*x1*x2*x3>>>