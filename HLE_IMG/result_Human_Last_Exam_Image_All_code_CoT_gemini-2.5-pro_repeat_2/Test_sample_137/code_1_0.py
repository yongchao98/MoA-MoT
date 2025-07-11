import numpy as np

def solve_ligand_problem():
    """
    Calculates the predicted distance between polymer chains for ligand 8.
    """
    # Step 1: Define the known data points
    # x: number of carbons in the R group
    # y: distance between adjacent polymer chains in Angstroms
    x_data = np.array([1, 2, 3])
    y_data = np.array([12.4, 12.6, 13.2])
    x_target = 8

    print("--- Problem Analysis ---")
    print(f"Data for Ligand 1 (C{x_data[0]}): {y_data[0]} Å")
    print(f"Data for Ligand 2 (C{x_data[1]}): {y_data[1]} Å")
    print(f"Data for Ligand 3 (C{x_data[2]}): {y_data[2]} Å")
    print(f"Goal: Predict distance for Ligand 8 (C{x_target}).\n")
    
    print("--- Modeling and Calculation ---")
    
    # Step 2: Model A - Linear extrapolation based on the overall trend (ligands 1 and 3)
    m_A = (y_data[2] - y_data[0]) / (x_data[2] - x_data[0])
    # Using point-slope form starting from ligand 1's data
    y_pred_A = y_data[0] + (x_target - x_data[0]) * m_A
    print("Model A (based on overall C1-C3 trend):")
    print(f"  Slope = ({y_data[2]} - {y_data[0]}) / ({x_data[2]} - {x_data[0]}) = {m_A:.2f} Å/carbon")
    print(f"  Prediction for C{x_target} = {y_data[0]} + ({x_target} - {x_data[0]}) * {m_A:.2f} = {y_pred_A:.1f} Å\n")

    # Step 3: Model B - Linear extrapolation based on the most recent trend (ligands 2 and 3)
    m_B = (y_data[2] - y_data[1]) / (x_data[2] - x_data[1])
    # Using point-slope form starting from ligand 3's data
    y_pred_B = y_data[2] + (x_target - x_data[2]) * m_B
    print("Model B (based on recent C2-C3 trend):")
    print(f"  Slope = ({y_data[2]} - {y_data[1]}) / ({x_data[2]} - {x_data[1]}) = {m_B:.2f} Å/carbon")
    print(f"  Prediction for C{x_target} = {y_data[2]} + ({x_target} - {x_data[2]}) * {m_B:.2f} = {y_pred_B:.1f} Å\n")
    
    # Step 4: Final prediction by averaging the results of the two models
    final_prediction = (y_pred_A + y_pred_B) / 2
    
    print("--- Final Result ---")
    print("Averaging the two plausible models gives a robust prediction.")
    print("Final Equation:")
    print(f"({y_pred_A:.1f} + {y_pred_B:.1f}) / 2 = {final_prediction:.1f}")
    
    print(f"\nThe predicted distance is {final_prediction:.1f} Å.")
    print("This corresponds to option C.")

solve_ligand_problem()