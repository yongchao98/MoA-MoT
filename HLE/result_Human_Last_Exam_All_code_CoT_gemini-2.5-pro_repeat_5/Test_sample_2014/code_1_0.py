from fractions import Fraction

def calculate_lift_ratio():
    """
    Calculates the lift ratio L1/L2 for two aerofoils in tandem
    in ground effect using the mirror image method.
    """
    print("Step 1: Define geometric parameters based on chord 'c'.")
    # We can set c=1 as the result is a non-dimensional ratio.
    # Using the Fraction class for exact arithmetic.
    c = Fraction(1)
    
    # Ride height h = c/2
    h = c / 2
    
    # Separation between aerofoils s = c/2
    s = c / 2
    
    # Horizontal distance between the vortices
    # d = (c/2) + s + (c/2) = c + s
    d = c + s
    
    print(f"Given chord c, ride height h = c/2, and separation s = c/2.")
    print(f"Let c = {c}.")
    print(f"This gives h = {h} and s = {s}.")
    print(f"The horizontal distance between vortices, d = c + s = {d}.")
    print("-" * 30)

    print("Step 2: Calculate the interaction parameter 'B'.")
    # The parameter B is derived from lifting-line theory.
    # B = (2 * c * h^2) / (d * (d^2 + 4 * h^2))
    
    h_sq = h**2
    d_sq = d**2
    
    numerator_B = 2 * c * h_sq
    denominator_B = d * (d_sq + 4 * h_sq)
    
    B = numerator_B / denominator_B
    
    print("The formula for B is: 2*c*h^2 / (d * (d^2 + 4*h^2))")
    print(f"Numerator = 2 * {c} * ({h})^2 = {numerator_B}")
    print(f"Denominator = {d} * (({d})^2 + 4 * ({h})^2) = {denominator_B}")
    print(f"B = {numerator_B} / {denominator_B} = {B}")
    print("-" * 30)

    print("Step 3: Calculate the lift ratio L1/L2.")
    # The lift ratio L1/L2 is given by (1 - B) / (1 + B)
    # assuming both aerofoils are at the same angle of attack.
    
    numerator_ratio = 1 - B
    denominator_ratio = 1 + B
    
    lift_ratio = numerator_ratio / denominator_ratio
    
    print("The formula for the lift ratio is: L1/L2 = (1 - B) / (1 + B)")
    print(f"The equation with the calculated value of B is:")
    print(f"L1/L2 = (1 - {B}) / (1 + {B})")
    print(f"L1/L2 = ({numerator_ratio}) / ({denominator_ratio})")
    print(f"L1/L2 = {lift_ratio}")
    print("-" * 30)
    
    print(f"The final lift ratio L1/L2 is {lift_ratio.numerator}/{lift_ratio.denominator}.")
    print(f"As a decimal, this is approximately {float(lift_ratio):.5f}.")

# Run the calculation
calculate_lift_ratio()

# The final answer in the required format
final_answer = Fraction(35, 43)
# print(f'<<<{final_answer.numerator}/{final_answer.denominator}>>>')