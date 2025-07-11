import numpy as np

def calculate_trace_of_covariance():
    """
    This function calculates the trace of the covariance matrix for the specified sampling procedure.
    The calculation is done analytically based on the statistical properties of the variables.

    The trace of the covariance of v is the same as the trace of the covariance of d
    because v is an orthogonal transformation (Householder reflection) of d.
    Trace(Cov(v)) = Trace(Cov(d))

    The vector d is constructed to be a unit vector, i.e., ||d||^2 = 1.
    This simplifies the trace calculation:
    Trace(Cov(d)) = E[||d||^2] - ||E[d]||^2 = 1 - ||E[d]||^2

    The expectation of d, E[d], has all components equal to zero except the first,
    due to the symmetric properties of the standard normal distribution used for vector c.
    So, ||E[d]||^2 = (E[d_1])^2.

    The first component d_1 is (a-b)/(a+b). With a ~ Gamma(alpha, theta) and b ~ Gamma(beta, theta),
    the variable z = a/(a+b) follows a Beta(alpha, beta) distribution.
    d_1 = a/(a+b) - b/(a+b) = z - (1-z) = 2z - 1.
    E[d_1] = 2 * E[z] - 1.
    The mean of a Beta(alpha, beta) distribution is alpha / (alpha + beta).

    Therefore, Trace(Cov(v)) = 1 - (2 * alpha / (alpha + beta) - 1)^2.
    """

    # Given parameters
    alpha = 3
    beta = 2
    
    # Calculate the expected value of d_1
    # E[Beta(alpha, beta)] = alpha / (alpha + beta)
    E_z = alpha / (alpha + beta)
    E_d1 = 2 * E_z - 1
    
    # Calculate the squared norm of the expected value of d
    # ||E[d]||^2 = (E[d_1])^2
    norm_E_d_squared = E_d1**2
    
    # Calculate the trace of the covariance matrix
    # Trace(Cov(v)) = 1 - ||E[d]||^2
    trace_cov_v = 1 - norm_E_d_squared
    
    # Print the calculation steps with numbers
    print(f"Given alpha = {alpha}, beta = {beta}")
    print(f"The expectation of the first component of d is E[d_1] = 2 * ({alpha}/({alpha}+{beta})) - 1 = {E_d1:.2f}")
    print(f"The squared L2 norm of the expectation of d is ||E[d]||^2 = ({E_d1:.2f})^2 = {norm_E_d_squared:.4f}")
    print("The trace of the covariance matrix of v is Trace(Cov(v)) = 1 - ||E[d]||^2")
    print(f"Trace = 1 - {norm_E_d_squared:.4f} = {trace_cov_v:.4f}")
    print(f"Final equation: {1} - ({E_d1})**2 = {trace_cov_v}")
    
    # As a fraction, the result is 1 - (1/5)^2 = 1 - 1/25 = 24/25 = 0.96
    
# Run the calculation
calculate_trace_of_covariance()
