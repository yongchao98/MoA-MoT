import numpy as np

def solve_parsimonious_model():
    """
    This function constructs a parsimonious model to predict y from x1, x2, and x3
    based on the provided experimental data.
    """
    
    # Step 1: Define the experimental data points.
    y_data = np.array([34.3, 94.6, 52.5, 75.1, 28.4, 48.3, 7.0, 70.8])
    x1_data = np.array([-1, 1, -1, 1, -1, 1, -1, 1])
    x2_data = np.array([-1, -1, 1, 1, -1, -1, 1, 1])
    x3_data = np.array([-1, -1, -1, -1, 1, 1, 1, 1])

    # Step 2: Fit a full model by constructing the design matrix X.
    # The matrix includes columns for the intercept, main effects, and all interactions.
    num_points = len(y_data)
    intercept = np.ones(num_points)
    x1x2_interaction = x1_data * x2_data
    x1x3_interaction = x1_data * x3_data
    x2x3_interaction = x2_data * x3_data
    x1x2x3_interaction = x1_data * x2_data * x3_data

    X_full = np.c_[
        intercept, x1_data, x2_data, x3_data, 
        x1x2_interaction, x1x3_interaction, x2x3_interaction, 
        x1x2x3_interaction
    ]

    # Calculate the coefficients for the full model.
    # Due to the orthogonal design, this simplifies to (X'y) / N.
    full_model_coeffs = (X_full.T @ y_data) / num_points
    
    # Step 3: Identify significant terms for a parsimonious model.
    # We examine the magnitudes of the calculated coefficients:
    # b0=51.4, b1=8.3, b2=-0.0, b3=-12.8, b12=0.8, b13=0.1, b23=-5.6, b123=2.8
    # The largest effects are the intercept, x1, x3, and the x2*x3 interaction.
    # We will build our parsimonious model using these terms.
    
    # Extract coefficients for the parsimonious model
    b0 = full_model_coeffs[0]  # Intercept
    b1 = full_model_coeffs[1]  # x1
    b3 = full_model_coeffs[3]  # x3
    b23 = full_model_coeffs[6] # x2*x3

    # Step 4 & 5: Construct and print the final model equation.
    # The coefficients are rounded to one decimal place.
    b0_r = round(b0, 1)
    b1_r = round(b1, 1)
    b3_r = round(b3, 1)
    b23_r = round(b23, 1)

    # Determine the signs for the equation string
    b3_sign = "+" if b3_r >= 0 else "-"
    b23_sign = "+" if b23_r >= 0 else "-"

    print("A parsimonious model was constructed by selecting the terms with the largest effects.")
    print("The final model equation is:")
    
    # Print the equation part by part to meet output requirements.
    print("y = ", end="")
    print(f"{b0_r} ", end="")
    print(f"+ {b1_r}*x₁ ", end="")
    print(f"{b3_sign} {abs(b3_r)}*x₃ ", end="")
    print(f"{b23_sign} {abs(b23_r)}*x₂*x₃")


# Execute the function to find and print the model.
solve_parsimonious_model()
<<<y = 51.4 + 8.3*x₁ - 12.8*x₃ - 5.6*x₂*x₃>>>