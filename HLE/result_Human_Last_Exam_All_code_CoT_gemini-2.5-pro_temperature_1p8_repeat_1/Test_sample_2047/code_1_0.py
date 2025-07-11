import math

# This script calculates the value of l(d, lambda) for given d and lambda.
# As per the problem description, d must be an integer >= 4 and lambda >= 1.
# We will use the example case of d=4 and lambda=1.

d_val = 4
lambda_val = 1.0

print(f"Calculating l(d, lambda) for d = {d_val}, lambda = {lambda_val}")
print("=" * 40)

# Validate inputs based on the function's domain
if not isinstance(d_val, int) or d_val < 4:
    print(f"Error: d must be an integer >= 4. Received d = {d_val}")
elif not isinstance(lambda_val, (int, float)) or lambda_val < 1:
    print(f"Error: lambda must be a number >= 1. Received lambda = {lambda_val}")
else:
    # Step 1: Calculate t1 and t2
    t1 = math.acos(math.sqrt(3.0 / d_val))
    t2 = math.acos(math.sqrt(2.0 / d_val))

    # Step 2: Calculate the two main terms of the equation for l(d, lambda)
    # Term 1 = (t2^2 - t1^2) / (2 * lambda)
    term1_numerator = t2**2 - t1**2
    term1 = term1_numerator / (2 * lambda_val)

    # Term 2 = (d - 2) * ln( (t2 * sin(t1)) / (t1 * sin(t2)) )
    sin_t1 = math.sin(t1)
    sin_t2 = math.sin(t2)
    
    log_argument = (t2 * sin_t1) / (t1 * sin_t2)
    term2 = (d_val - 2) * math.log(log_argument)
    
    # Step 3: Sum the terms to get the final result
    final_ell_value = term1 + term2

    # Step 4: Print the full calculation breakdown
    print("The equation is: l(d, lambda) = (t2^2 - t1^2)/(2*lambda) + (d-2)*ln((t2*sin(t1))/(t1*sin(t2)))")
    
    print("\n--- Intermediate Values ---")
    print(f"t1 = arccos(sqrt(3/{d_val})) = arccos({math.sqrt(3.0 / d_val):.4f}) = {t1:.4f} radians")
    # For d=4, t1 = pi/6
    print(f"t2 = arccos(sqrt(2/{d_val})) = arccos({math.sqrt(2.0 / d_val):.4f}) = {t2:.4f} radians")
    # For d=4, t2 = pi/4

    print("\n--- Equation Term by Term ---")
    print(f"Term 1 = (t2^2 - t1^2) / (2 * lambda) = ({t2**2:.4f} - {t1**2:.4f}) / (2 * {lambda_val}) = {term1:.4f}")
    
    print(f"Term 2 = (d-2)*ln(...) = ({d_val-2}) * ln(({t2:.4f} * {sin_t1:.4f}) / ({t1:.4f} * {sin_t2:.4f})) = {term2:.4f}")

    print("\n--- Final Equation ---")
    print(f"l({d_val}, {lambda_val}) = {term1:.4f} + {term2:.4f}")
    print(f"l({d_val}, {lambda_val}) = {final_ell_value:.4f}")

    # For the case d=4, lambda=1:
    # t1 = pi/6, t2 = pi/4
    # t1^2 = pi^2/36, t2^2 = pi^2/16
    # sin(t1) = 1/2, sin(t2) = 1/sqrt(2)
    # term1 = (pi^2/16 - pi^2/36)/(2) = 5*pi^2/288
    # term2 = 2 * ln( (pi/4 * 1/2) / (pi/6 * 1/sqrt(2)) ) = 2 * ln( (6*sqrt(2))/8 ) = 2 * ln(3*sqrt(2)/4) = ln(18/16) = ln(9/8)
    exact_value_str = f"(5 * pi^2 / 288) + ln(9/8) = {(5*math.pi**2/288) + math.log(9/8):.10f}"
    print(f"\nExact value for d=4, lambda=1 is: {exact_value_str}")
