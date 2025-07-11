def print_magnetic_field_solution():
    """
    This function prints the derived symbolic expression for the magnetic field
    H(r, theta) inside and outside the spherical shell, corresponding to the correct answer choice.
    """
    # Unicode characters for mathematical symbols
    mu = "\u03BC"
    mu_0 = "\u03BC\u2080"
    theta = "\u03B8"
    phi_hat = "\u03C6\u0302"
    r_hat = "r\u0302"
    theta_hat = "\u03B8\u0302"
    z_hat = "z\u0302"
    infinity = "\u221e"
    
    print("The final expression for the magnetic field H(r, \u03B8) is:")
    print(f"H(r, {theta}) =")
    print("{")
    
    # Inside the sphere (0 < r < R)
    inside_field = f"    (2 * {mu_0} / {mu}) * (K\u2080 / (1 + 2 * {mu_0} / {mu})) {z_hat}, for 0 < r < R"
    print(inside_field)
    
    # Outside the sphere (R < r < infinity)
    outside_field_part1 = f"    (K\u2080 / (1 + 2 * {mu_0} / {mu}))"
    outside_field_part2 = f" * (R\u00b3 / r\u00b3) * (2*cos({theta})*{r_hat} + sin({theta})*{theta_hat}), for R < r < {infinity}"
    print(outside_field_part1 + outside_field_part2)
    
    print("}")

print_magnetic_field_solution()