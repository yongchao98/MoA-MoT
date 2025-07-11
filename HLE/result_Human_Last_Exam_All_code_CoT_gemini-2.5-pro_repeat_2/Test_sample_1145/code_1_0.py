import math

def calculate_particle_in_box_probability():
    """
    Calculates the probability of finding a particle in a specific interval
    of a 1D box for a given quantum state n.
    """
    # 1. Define the constants for the problem
    n = 2
    x1_norm = 0.495  # Lower bound of the interval, x1/a
    x2_norm = 0.505  # Upper bound of the interval, x2/a

    # 2. Define a function for the antiderivative of the probability density.
    # The integral of (2/a)*sin^2(n*pi*x/a) dx is x/a - sin(2*n*pi*x/a)/(2*n*pi)
    def antiderivative(x_norm, n_val):
        """Calculates the value of the antiderivative F(x/a)."""
        constant = 2 * n_val * math.pi
        angle = constant * x_norm
        return x_norm - (1 / constant) * math.sin(angle)

    # 3. Calculate the value of the antiderivative at the interval boundaries
    val_at_x2 = antiderivative(x2_norm, n)
    val_at_x1 = antiderivative(x1_norm, n)

    # 4. The probability is the difference between the values at the boundaries
    probability = val_at_x2 - val_at_x1
    
    # 5. Output the final equation with all numbers and the result
    print("The probability (P) is the definite integral of the probability density function.")
    print(f"For n={n}, the integral is evaluated from x={x1_norm}a to x={x2_norm}a.")
    print("\nThe final equation is:")
    
    # Constructing the string for the equation as requested
    term1_n = 2 * n
    term1_x = x2_norm
    term2_n = 2 * n
    term2_x = x1_norm
    
    equation_str = (f"P = ({term1_x} - sin({term1_n}*pi*{term1_x}) / ({term1_n}*pi)) - "
                    f"({term2_x} - sin({term2_n}*pi*{term2_x}) / ({term2_n}*pi))")
    
    print(equation_str)
    
    # Print the calculated result
    print(f"\nWhich evaluates to:")
    print(f"P = ({val_at_x2}) - ({val_at_x1})")
    print(f"\nFinal Probability = {probability}")

# Execute the function
calculate_particle_in_box_probability()
<<<6.579543163351322e-06>>>