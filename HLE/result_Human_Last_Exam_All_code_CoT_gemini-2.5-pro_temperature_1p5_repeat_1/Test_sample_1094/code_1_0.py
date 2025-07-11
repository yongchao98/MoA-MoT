import math

def calculate_normalized_loss(i):
    """
    Calculates the normalized AC loss for an elliptical superconductor.
    The formula is 2 * [(1-i)ln(1-i) + (1+i)ln(1+i) - i^2].
    """
    # The limit of x*ln(x) as x->0 is 0.
    # Handle the terms separately to avoid issues at the boundaries if necessary.
    # However, since the prompt is for i < 1, direct calculation is fine.
    
    term1 = (1 - i) * math.log(1 - i)
    term2 = (1 + i) * math.log(1 + i)
    term3 = i**2
    
    result = 2 * (term1 + term2 - term3)
    return result

def main():
    """
    Main function to get user input and display the results.
    """
    i = -1.0
    while True:
        try:
            i_str = input("Enter the value of the normalized current i = Im/Ic (must be between 0 and 1): ")
            i = float(i_str)
            if 0 < i < 1:
                break
            else:
                print("Error: The value of i must be strictly between 0 and 1.")
        except ValueError:
            print("Invalid input. Please enter a number.")
        except Exception as e:
            print(f"An error occurred: {e}")
            return
            
    # Calculate the normalized loss
    normalized_loss = calculate_normalized_loss(i)
    
    # Output the formula and the calculation
    print("\nThe normalized AC loss is given by the formula: 2πQ / (μ₀ * Ic²) = 2 * [(1-i)ln(1-i) + (1+i)ln(1+i) - i²]")
    print("\nFor your input value i =", i)
    
    # Print the equation with the numbers substituted
    print("\nCalculation:")
    val_1_minus_i = 1 - i
    val_1_plus_i = 1 + i
    val_i_squared = i**2
    
    # We use .4f formatting for better readability of intermediate steps
    # but the actual calculation uses full precision
    print(f"2 * [(1 - {i}) * ln(1 - {i}) + (1 + {i}) * ln(1 + {i}) - {i}²]")
    print(f"= 2 * [({val_1_minus_i:.4f}) * ln({val_1_minus_i:.4f}) + ({val_1_plus_i:.4f}) * ln({val_1_plus_i:.4f}) - {val_i_squared:.4f}]")
    
    # Print the final result
    print(f"\nFinal Result for 2πQ / (μ₀ * Ic²):")
    print(normalized_loss)
    
    # Final answer in the specified format
    print(f"<<<{normalized_loss}>>>")

if __name__ == "__main__":
    main()