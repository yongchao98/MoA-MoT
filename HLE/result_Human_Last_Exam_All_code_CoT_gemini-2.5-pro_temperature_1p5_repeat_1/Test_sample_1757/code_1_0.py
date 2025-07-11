import math

def calculate_sum():
    """
    Calculates the sum over all natural dimensions of the specified ratio.
    """
    
    # Calculate the constant C_dist = E[sqrt(t_i^2 + t_j^2)], where t_i, t_j ~ U(-1, 1).
    # The analytical result is (sqrt(2) + ln(1 + sqrt(2))) / 3.
    C_dist = (math.sqrt(2) + math.log(1 + math.sqrt(2))) / 3.0

    total_sum = 0.0
    
    print("Calculating the sum of the series S = R_1 + R_2 + R_3 + ...")
    print("where R_n = (n+1) / (n! * 2^n * (1 + (n-1)*C_dist))\n")
    
    equation_str = "S = "
    
    # The series converges rapidly. A loop up to n=20 is more than sufficient
    # for precision well beyond the third decimal place.
    for n in range(1, 21):
        try:
            # Numerator of the n-th ratio R_n
            numerator_Rn = n + 1
            
            # Denominator of the n-th ratio R_n
            denominator_Rn = (math.factorial(n) * (2**n) * 
                            (1 + (n - 1) * C_dist))
            
            term = numerator_Rn / denominator_Rn
        except OverflowError:
            # The term becomes too small to calculate with standard floats,
            # so we can stop summing.
            break
            
        # Add the numeric value of the term to the equation string
        if n > 1:
            equation_str += " + "
        equation_str += f"{term:.8f}"

        total_sum += term
        
        # Stop if the term is too small to affect the desired precision.
        if term < 1e-12:
            break
    
    print("The final equation with its calculated numbers is:")
    print(equation_str + " + ...")
    
    print(f"\nFinal sum: {total_sum}")
    print(f"The sum determined with three-decimal precision is: {total_sum:.3f}")

# Execute the calculation
calculate_sum()