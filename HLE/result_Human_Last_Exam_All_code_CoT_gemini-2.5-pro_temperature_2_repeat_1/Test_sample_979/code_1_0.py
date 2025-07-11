def print_magnetic_field_solution():
    """
    Prints the solution for the magnetic field H inside and outside the sphere.
    """
    
    # These expressions correspond to option E, which is the physically correct result,
    # likely assuming a sign convention difference in the problem statement.
    
    h_in_expression = "H_in(r, theta) = (2 * mu_0 * K_0) / (mu * (1 + (2 * mu_0) / mu)) * z_hat"
    h_out_expression = "H_out(r, theta) = (K_0 * R^3) / (r^3 * (1 + (2 * mu_0) / mu)) * (2 * cos(theta) * r_hat + sin(theta) * theta_hat)"

    print("The magnetic field H(r, theta) is given by:")
    print("For 0 < r < R:")
    print(f"  {h_in_expression}")
    print("\nFor r > R:")
    print(f"  {h_out_expression}")
    
    print("\nWhich can be written as:")
    
    final_answer = r"""
H(r, \theta) = 
    \begin{cases} 
    \frac{2 \mu_0}{\mu} \frac{K_0}{1 + \frac{2 \mu_0}{\mu}} \hat{z} & 0 < r < R \\
    \frac{K_0}{1 + \frac{2 \mu_0}{\mu}} \frac{R^3}{r^3} (2 \cos \theta \hat{r} + \sin \theta \hat{\theta}) & R < r < \infty
    \end{cases}
"""
    print(final_answer)

print_magnetic_field_solution()