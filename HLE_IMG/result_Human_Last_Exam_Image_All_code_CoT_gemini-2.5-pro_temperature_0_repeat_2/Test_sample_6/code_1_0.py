def solve_histopathology_task():
    """
    Analyzes statements about kidney biopsy images to identify the correct ones.
    This function codifies the reasoning based on visual analysis of the provided images.
    """
    
    # A detailed analysis determines the truth value of each statement.
    # Statement 1: False. Image D shows prominent Kimmelstiel-Wilson lesions (nodular glomerulosclerosis).
    # Statement 2: True. Image C shows classic features of nodular glomerulosclerosis.
    # Statement 3: False. The described feature, "effacement of Bowman's capsule," is not clearly visible in Image C.
    # Statement 4: True. The arrows in Image D correctly point to the extracellular matrix deposits forming the nodules.
    
    true_statements = [2, 4]
    
    print("Step-by-step analysis of the statements:")
    print("1. Statement 1 is evaluated as FALSE. While Image B shows segmental sclerosis, Image D exhibits prominent nodular glomerulosclerosis (Kimmelstiel-Wilson lesions), contradicting the claim that they are absent.")
    print("2. Statement 2 is evaluated as TRUE. Image C is a classic example of eosinophilic nodular glomerulosclerosis with mesangial matrix expansion.")
    print("3. Statement 3 is evaluated as FALSE. 'Effacement of Bowman's capsule' is an imprecise term for light microscopy, and this feature is not clearly observed in Image C, where Bowman's space remains patent.")
    print("4. Statement 4 is evaluated as TRUE. The arrows in Image D correctly indicate the large deposits of extracellular matrix that are the hallmark of nodular glomerulosclerosis.")
    
    print("\nConclusion:")
    print("The statements that are true are:")
    # Printing each number of the true statements as requested.
    for number in true_statements:
        print(f"Statement {number}")
        
    # The combination [2, 4] corresponds to answer choice G.
    print("\nThis set of true statements, [2, 4], corresponds to answer choice G.")

solve_histopathology_task()