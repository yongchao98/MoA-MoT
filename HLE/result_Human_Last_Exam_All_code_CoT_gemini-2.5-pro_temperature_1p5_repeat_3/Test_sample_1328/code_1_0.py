import numpy as np

def solve_optimal_feedback():
    """
    This function calculates the optimal feedback adjustment factor 'b' for a
    wireless communication system.
    
    The user is prompted to enter the required parameters.
    """
    try:
        # User-provided values for alpha and P1.
        # Example values can be alpha = 0.5, P1 = 1.0
        alpha_str = input("Enter the weather-induced correlation (alpha, e.g., 0.5): ")
        alpha = float(alpha_str)
        if not (-1 < alpha < 1):
            print("Error: alpha must be between -1 and 1 for a valid covariance matrix.")
            return

        p1_str = input("Enter the power used in the first transmission (P1, e.g., 1.0): ")
        p1 = float(p1_str)
        if p1 < 0:
            print("Error: Power P1 cannot be negative.")
            return

        # The optimal feedback factor b is calculated using the derived formula.
        # b = -alpha * P1
        b_opt = -alpha * p1

        # Output the derivation and the final result
        print("\n--- Derivation ---")
        print("The mutual information is maximized by maximizing the determinant |K_X+Z|.")
        print("The part of the determinant depending on 'b' is f(b) = -b^2 - 2*alpha*P1*b.")
        print("This is a downward parabola in 'b', with its maximum at the vertex.")
        print("The vertex is found by setting the derivative f'(b) = -2b - 2*alpha*P1 to zero.")
        print("Solving for b gives the optimal feedback factor.")
        
        print("\n--- Optimal Feedback Factor ---")
        print(f"Given alpha = {alpha} and P1 = {p1}:")
        print("The optimal feedback factor b is calculated as:")
        # To display the full equation, we handle the sign of alpha
        if alpha >= 0:
            print(f"b = - ( {alpha} * {p1} )")
        else:
            print(f"b = - ( ({alpha}) * {p1} )")

        print(f"b = {b_opt}")
        
        print("\nNote: This solution assumes it is feasible within the power constraints of the second transmission (i.e., b^2 <= P2).")

    except ValueError:
        print("Invalid input. Please enter numerical values.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Run the function
solve_optimal_feedback()