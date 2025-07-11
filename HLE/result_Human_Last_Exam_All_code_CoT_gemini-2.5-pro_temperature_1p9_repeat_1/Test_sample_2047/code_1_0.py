import math

def calculate_l(d, lambda_val):
    """
    Calculates the value of l(d, lambda) based on the derived formula.
    
    The formula is derived under the assumption that p_{d,λ}(x) is proportional to
    exp(-θ²/(2λ)) * (sin(θ)/θ)^(d-2), where θ is the angular distance of x
    from a pole p₀ = 1_d/√d.
    
    l(d,λ) = ln[ p(x₁) / p(x₂) ]
           = (θ₂² - θ₁²)/(2λ) + (d-2) * ln( (sin(θ₁)/θ₁) / (sin(θ₂)/θ₂) )
           
    where:
    p₀ = (1/√d, ..., 1/√d)
    x₁ = (e₁ + e₂ + e₃)/√3  => x₁⋅p₀ = √3/√d
    x₂ = (e₃ + e₄)/√2      => x₂⋅p₀ = √2/√d
    θ₁ = arccos(√(3/d))
    θ₂ = arccos(√(2/d))
    """
    
    # Check for domain validity
    if d < 4:
        print("Error: d must be >= 4")
        return None
    if lambda_val < 1:
        print("Error: lambda must be >= 1")
        return None

    # Calculate cos(theta) values
    cos_theta1_sq = 3 / d
    cos_theta2_sq = 2 / d
    
    if cos_theta1_sq > 1 or cos_theta2_sq > 1:
        print("Error: d is too small, resulting in arccos argument > 1.")
        return None

    cos_theta1 = math.sqrt(cos_theta1_sq)
    cos_theta2 = math.sqrt(cos_theta2_sq)
    
    # Calculate theta values
    theta1 = math.acos(cos_theta1)
    theta2 = math.acos(cos_theta2)
    
    # Calculate the first term: (θ₂² - θ₁²) / (2λ)
    term1 = (theta2**2 - theta1**2) / (2 * lambda_val)
    
    # Calculate the second term: (d-2) * ln( (sin(θ₁)/θ₁) / (sin(θ₂)/θ₂) )
    # Handle theta=0 case to avoid division by zero, although not possible for d>=4.
    sinc1 = math.sin(theta1) / theta1 if theta1 != 0 else 1.0
    sinc2 = math.sin(theta2) / theta2 if theta2 != 0 else 1.0
    
    term2 = (d - 2) * math.log(sinc1 / sinc2)
    
    # Final result
    result = term1 + term2
    
    # Print the equation and intermediate values as requested
    print("Calculating l(d, λ) for d = {} and λ = {}".format(d, lambda_val))
    print("-" * 50)
    print("Formula: l(d, λ) = (θ₂² - θ₁²)/(2λ) + (d-2) * ln(sinc(θ₁)/sinc(θ₂))")
    print("where θ₁ = arccos(√({}/{})) and θ₂ = arccos(√({}/{}))".format(3,d,2,d))
    print("-" * 50)
    
    print(f"θ₁ = arccos({cos_theta1:.8f}) = {theta1:.8f} radians")
    print(f"θ₂ = arccos({cos_theta2:.8f}) = {theta2:.8f} radians")
    print("\n--- Term 1 Calculation ---")
    print(f"(θ₂² - θ₁²) / (2λ) = (({theta2:.8f})² - ({theta1:.8f})²) / (2 * {lambda_val})")
    print(f"                 = {theta2**2 - theta1**2:.8f} / {2*lambda_val:.1f}")
    print(f"                 = {term1:.8f}")
    
    print("\n--- Term 2 Calculation ---")
    print(f"sinc(θ₁) = sin({theta1:.8f}) / {theta1:.8f} = {sinc1:.8f}")
    print(f"sinc(θ₂) = sin({theta2:.8f}) / {theta2:.8f} = {sinc2:.8f}")
    print(f"(d-2) * ln(sinc(θ₁)/sinc(θ₂)) = ({d}-2) * ln({sinc1:.8f}/{sinc2:.8f})")
    print(f"                             = {d-2} * ln({sinc1/sinc2:.8f})")
    print(f"                             = {d-2} * {math.log(sinc1/sinc2):.8f}")
    print(f"                             = {term2:.8f}")
    
    print("\n--- Final Result ---")
    print(f"l({d}, {lambda_val}) = {term1:.8f} + {term2:.8f}")
    print(f"l({d}, {lambda_val}) = {result:.8f}")
    
    return result

if __name__ == '__main__':
    # Using example values d=4, λ=1 as they are the smallest integers in the function's domain.
    d_val = 4
    lambda_val = 1
    final_answer = calculate_l(d_val, lambda_val)
