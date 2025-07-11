import math

def calculate_lower_bound():
    """
    This function calculates the constant lower bound for d(t,x).
    The lower bound L is determined by the minimum value of a function L(u)
    for u in the interval [0, 1].
    The analytical derivation shows that the minimum is attained at u = 1.
    """
    
    # The value of u where the minimum of the lower bound function is achieved
    u = 1.0

    print("The constant lower bound L is calculated from the formula:")
    print("L = ( (3*u - 5*u^2) - sqrt((3*u - 5*u^2)^2 + 8*u^3*(1-u)) ) / 4")
    print(f"We evaluate this formula at u = {u}, where the minimum is attained.")
    print("\n--- Calculation Steps ---")
    
    # Numerator part 1: (3*u - 5*u^2)
    term1_val = 3 * u - 5 * u**2
    print(f"Step 1: Calculate the first term in the numerator: (3*{u} - 5*{u}^2) = {term1_val}")
    
    # Numerator part 2: sqrt(...)
    # Base of the square root
    sqrt_base_term1_val = (3 * u - 5 * u**2)**2
    sqrt_base_term2_val = 8 * u**3 * (1 - u)
    sqrt_base_val = sqrt_base_term1_val + sqrt_base_term2_val
    print(f"Step 2: Calculate the term under the square root:")
    print(f"  ({term1_val})^2 + 8*({u})^3*(1-{u}) = {sqrt_base_term1_val} + {sqrt_base_term2_val} = {sqrt_base_val}")
    
    # Square root value
    sqrt_val = math.sqrt(sqrt_base_val)
    print(f"Step 3: Calculate the square root: sqrt({sqrt_base_val}) = {sqrt_val}")
    
    # Full numerator
    numerator_val = term1_val - sqrt_val
    print(f"Step 4: Calculate the full numerator: {term1_val} - {sqrt_val} = {numerator_val}")
    
    # Denominator
    denominator_val = 4.0
    
    # Final result
    L = numerator_val / denominator_val
    print(f"Step 5: Divide by the denominator: {numerator_val} / {denominator_val} = {L}")
    
    print("\n--- Final Equation with Numbers ---")
    print(f"L = (({term1_val}) - sqrt(({term1_val})^2 + ({sqrt_base_term2_val}))) / {denominator_val}")
    print(f"L = (({term1_val}) - sqrt({sqrt_base_term1_val} + {sqrt_base_term2_val})) / {denominator_val}")
    print(f"L = (({term1_val}) - sqrt({sqrt_base_val})) / {denominator_val}")
    print(f"L = (({term1_val}) - {sqrt_val}) / {denominator_val}")
    print(f"L = ({numerator_val}) / {denominator_val}")
    print(f"L = {L}")
    
    return L

if __name__ == '__main__':
    lower_bound = calculate_lower_bound()
    print(f"\nThe constant lower bound of d(t,x) is: {lower_bound}")
    print(f"<<<{-1.0}>>>")