from fractions import Fraction

def solve_conductance_ratio():
    """
    Calculates the ratio of the fourth to the first moment of conductance
    for a critical disordered Majorana wire (class D).
    """
    
    # Step 1: Define the known seed values for the moments of the dimensionless conductance g.
    # The average conductance <g> for class D is a universal value.
    g1 = Fraction(1, 2)
    
    # The variance var(g) = <g^2> - <g>^2 is also universal.
    var_g = Fraction(1, 16)
    
    # From this, we can find the second moment <g^2>.
    g2 = var_g + g1**2
    
    # Step 2: Use the exact recurrence relation for the moments of g in class D
    # for n >= 2:  2(n-1)<g^(n+1)> = (n-2)<g^n> + <g^(n-1)>
    
    # Calculate <g^3> by setting n=2 in the recurrence relation:
    # 2(2-1)<g^3> = (2-2)<g^2> + <g^1>
    # 2<g^3> = <g^1>
    g3 = g1 / 2
    
    # Calculate <g^4> by setting n=3 in the recurrence relation:
    # 2(3-1)<g^4> = (3-2)<g^3> + <g^2>
    # 4<g^4> = <g^3> + <g^2>
    g4 = (g3 + g2) / 4
    
    # Step 3: Calculate the required ratio of the fourth moment to the average value.
    ratio = g4 / g1
    
    # Step 4: Print the final equation with each number, as requested.
    print("For a critical disordered Majorana wire:")
    print(f"The average conductance is <g> = {g1.numerator}/{g1.denominator}")
    print(f"The fourth statistical moment of conductance is <g^4> = {g4.numerator}/{g4.denominator}")
    print("\nThe ratio <g^4> / <g> is calculated as:")
    print(f"({g4.numerator}/{g4.denominator}) / ({g1.numerator}/{g1.denominator}) = {ratio.numerator}/{ratio.denominator}")
    
    # Print the final numerical answer
    print(f"\nThe value of the ratio is: {float(ratio)}")

solve_conductance_ratio()