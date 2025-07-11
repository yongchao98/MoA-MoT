from fractions import Fraction

def calculate_variance():
    """
    Calculates the variance of Y based on the analytical solution
    using order statistics.
    """
    n = 4

    # E[Y] calculation
    # E[Y] = 1/4 * E[X(3)] +
    #        1/4 * (1/2 * E[X(1)] + 1/2 * E[X(3)]) +
    #        1/4 * (1/2 * E[X(2)] + 1/2 * E[X(4)]) +
    #        1/4 * E[X(2)]
    
    # E[X(j)] = j / (n+1)
    e_x = [Fraction(j, n + 1) for j in range(n + 2)] # e_x[j] = E[X(j)]
    
    e_y = Fraction(1, 4) * e_x[3] + \
          Fraction(1, 8) * (e_x[1] + e_x[3]) + \
          Fraction(1, 8) * (e_x[2] + e_x[4]) + \
          Fraction(1, 4) * e_x[2]
          
    # E[Y^2] calculation
    # E[Y^2] = 1/4 * E[X(3)^2] +
    #          1/4 * (1/2 * E[X(1)^2] + 1/2 * E[X(3)^2]) +
    #          1/4 * (1/2 * E[X(2)^2] + 1/2 * E[X(4)^2]) +
    #          1/4 * E[X(2)^2]
    # This simplifies to:
    # E[Y^2] = 1/8*E[X(1)^2] + 3/8*E[X(2)^2] + 3/8*E[X(3)^2] + 1/8*E[X(4)^2]

    # E[X(j)^2] = j*(j+1) / ((n+1)*(n+2))
    e_x_sq = [Fraction(j * (j + 1), (n + 1) * (n + 2)) for j in range(n + 2)] # e_x_sq[j] = E[X(j)^2]
    
    e_y_sq = Fraction(1, 8) * e_x_sq[1] + \
             Fraction(3, 8) * e_x_sq[2] + \
             Fraction(3, 8) * e_x_sq[3] + \
             Fraction(1, 8) * e_x_sq[4]

    # Variance Calculation: Var(Y) = E[Y^2] - (E[Y])^2
    var_y = e_y_sq - e_y**2
    
    print("The expectation of Y is E[Y] = {}".format(e_y))
    print("The expectation of Y^2 is E[Y^2] = {}".format(e_y_sq))
    print("\nThe variance of Y is given by the equation:")
    print("Var(Y) = E[Y^2] - (E[Y])^2")
    print("Var(Y) = {} - ({})^2 = {} - {}".format(e_y_sq, e_y, e_y_sq, e_y**2))
    print("Var(Y) = {}".format(var_y))

calculate_variance()