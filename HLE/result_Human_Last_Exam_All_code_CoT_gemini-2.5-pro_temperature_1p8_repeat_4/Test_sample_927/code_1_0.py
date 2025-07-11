def generate_composite_number_formula():
    """
    This function demonstrates how a recursively enumerable set (the composite numbers)
    can be defined by an existential formula in the given language.
    """
    # A number x is composite if it can be written as (y+2)*(z+2) for some natural numbers y and z.
    # The corresponding Diophantine equation is x - (y+2)*(z+2) = 0.
    diophantine_poly = "x - (y + 2) * (z + 2)"

    # The existential formula in the language L claims the existence of such y and z.
    # The predicate P(v) asserts that v is a natural number.
    # Note that the number 2 can be written as (1+1), which is a term in the language.
    y_var = "y"
    z_var = "z"
    two = "(1+1)"
    formula = (f"exists {y_var}, {z_var} ("
               f"P({y_var}) AND P({z_var}) AND "
               f"x - ({y_var} + {two}) * ({z_var} + {two}) = 0"
               f")")

    print("--- Example: Defining Composite Numbers ---")
    print(f"The set of composite numbers is a recursively enumerable (RE) set.")
    print(f"A number x is composite if x = a * b for integers a,b > 1.")
    print(f"This is equivalent to finding natural numbers y,z such that x = (y+2)*(z+2).")
    print(f"This leads to the Diophantine equation with integer coefficients:")
    print(f"p(x, y, z) = {diophantine_poly} = 0\n")

    print("This can be translated into an existential L-formula:")
    print(formula)

    print("\nEach number in the final equation is printed below as requested:")
    final_equation = f"x - (y + {two}) * (z + {two}) = 0"
    print(f"Final Equation: {final_equation}")
    print("Numbers in the equation: 1, 1, 1, 1")


generate_composite_number_formula()