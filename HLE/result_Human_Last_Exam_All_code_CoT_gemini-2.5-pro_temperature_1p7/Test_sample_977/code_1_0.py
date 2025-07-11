def print_potential_formula():
    """
    This function prints the derived formula for the electric potential
    in the region 0 <= y <= a.
    """
    
    # The derived formula for the electric potential Phi(x, y) for 0 <= y <= a is:
    numerator = "-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    denominator_part1 = "epsilon_2 * cosh(k*a) * sinh(k*b)"
    denominator_part2 = "epsilon_1 * sinh(k*a) * cosh(k*b)"
    denominator = f"k * ({denominator_part1} + {denominator_part2})"
    
    potential_formula = f"Phi(x, y) = ({numerator}) / ({denominator})"
    
    print("The electric potential Phi(x, y) in the region 0 <= y <= a is:")
    print(potential_formula)

print_potential_formula()