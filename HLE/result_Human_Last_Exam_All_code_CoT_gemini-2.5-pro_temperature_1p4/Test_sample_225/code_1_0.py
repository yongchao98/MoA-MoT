def solve_lissajous_intersections():
    """
    Calculates the number of self-intersection points for the curve
    x(t) = cos(9t), y(t) = sin(5t).
    """
    # The parameters p and q are the frequencies from the parametric equations.
    p = 9
    q = 5

    # The formula for the number of self-intersections (N) for a Lissajous curve
    # of the form x=cos(pt), y=sin(qt) where p and q are odd coprime integers is:
    # N = (p - 1) * (q - 1) / 2
    
    # We will now calculate the result and print each step as requested.
    
    p_minus_1 = p - 1
    q_minus_1 = q - 1
    product = p_minus_1 * q_minus_1
    num_intersections = product / 2

    print("The number of self-intersection points is calculated using the formula N = (p - 1) * (q - 1) / 2.")
    print("For p = 9 and q = 5, the equation is:")
    # The request is to output each number in the final equation.
    print(f"N = ({p} - 1) * ({q} - 1) / 2")
    print(f"N = {p_minus_1} * {q_minus_1} / 2")
    print(f"N = {product} / 2")
    print(f"N = {int(num_intersections)}")
    print("\nThe total number of self-intersection points is:")
    print(int(num_intersections))

solve_lissajous_intersections()