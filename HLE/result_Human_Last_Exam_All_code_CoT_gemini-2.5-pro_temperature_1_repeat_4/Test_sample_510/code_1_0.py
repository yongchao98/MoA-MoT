import sympy

def solve_queueing_problem():
    """
    Calculates the optimal mean response time for the given M/G/1 SRPT queue
    and identifies the components of the final expression.
    """
    # 1. Define the parameters of the queueing system
    # Arrival rate lambda
    lambda_val = sympy.Rational(3, 2)
    # Job size distribution is Uniform(0, 1), so f_S(s) = 1 for s in [0, 1]
    f_s_dist = 1
    # The integration range for job sizes is from 0 to 1
    s_min, s_max = 0, 1

    # Define symbolic variables for integration
    s, y = sympy.symbols('s y')

    # 2. Define the formula for the optimal mean response time x (SRPT policy)
    # x = integral_{s_min}^{s_max} (s + E[W_s]) * f_S(s) ds
    # where E[W_s] is the mean waiting time for a job of size s.
    # E[W_s] = (lambda * integral_0^s y^2*f_S(y)dy) / (2 * (1 - rho_s)^2)
    # and rho_s = lambda * integral_0^s y*f_S(y)dy

    # 3. Calculate the components of the formula
    # rho_s: traffic intensity from jobs of size up to s
    rho_s = lambda_val * sympy.integrate(y * f_s_dist, (y, 0, s))

    # Numerator of E[W_s]
    wait_time_numerator = lambda_val * sympy.integrate(y**2 * f_s_dist, (y, 0, s))

    # Denominator of E[W_s]
    wait_time_denominator = 2 * (1 - rho_s)**2

    # Mean waiting time for a job of size s
    E_W_s = wait_time_numerator / wait_time_denominator

    # 4. Calculate the two parts of the total mean response time integral
    
    # First part: Mean service time E[S] = integral(s * f_S(s) ds)
    mean_service_time = sympy.integrate(s * f_s_dist, (s, s_min, s_max))
    
    # Second part: Mean waiting time E[W] = integral(E[W_s] * f_S(s) ds)
    mean_waiting_time = sympy.integrate(E_W_s * f_s_dist, (s, s_min, s_max))

    # The optimal mean response time x is the sum
    x = mean_service_time + mean_waiting_time
    
    # 5. Output the results and the numbers in the final equation
    print("Step-by-step calculation of the optimal mean response time x:")
    print(f"Mean service time E[S] = {mean_service_time}")
    print(f"Mean waiting time E[W] = {mean_waiting_time}")
    print(f"Total optimal mean response time x = E[S] + E[W] = {x}")
    print("\nThe final expression for x is of the form: R + C * ln(A)")
    
    # Extract the components of the expression for x
    # x.as_coeff_add() separates the expression into a constant term and the rest
    rational_part, remaining_part = x.as_coeff_add()
    
    # remaining_part.as_coeff_Mul() separates a term into its coefficient and the main function
    coeff, log_term = remaining_part.as_coeff_Mul()
    
    # log_term.args[0] gives the argument of the logarithm
    log_argument = log_term.args[0]
    
    print("We identify the components of this expression:")
    print(f"The additive rational term (R) is: {rational_part}")
    print(f"The coefficient of the logarithm (C) is: {coeff}")
    print(f"The argument of the logarithm (A) is: {log_argument}")
    
    print("\nThe problem asks for the term remaining after removing rational terms and logarithms of rational numbers.")
    print(f"The rational term to remove is {rational_part}.")
    print(f"The remaining term is {remaining_part}.")
    print(f"This term can be written as ln({log_argument}**({coeff})).")
    check_term = log_argument**coeff
    # A number is rational in sympy if its is_rational property is True
    is_log_of_rational = check_term.is_rational
    print(f"The base of this logarithm, {check_term}, is rational: {is_log_of_rational}.")
    if not is_log_of_rational:
      print("Therefore, this term is not a logarithm of a rational number and should not be removed.")
      print(f"The final remaining term is: {remaining_part}")

if __name__ == '__main__':
    solve_queueing_problem()