def find_missing_boxes():
    original = [87, 123, 302, 308, 370, 86, 182, 79, 135, 170, 198, 243, 45, 212, 183, 174, 335, 134, 58, 131, 124, 360, 192, 245]
    solution = [
        370, 308, 360, 302, 335, 245, 58, 243, 212, 135, 198, 192, 134,
        183, 182, 131, 174, 170, 124, 45, 123, 87, 86, 79
    ]
    
    original_set = set(original)
    solution_set = set(solution)
    
    missing = original_set - solution_set
    print("Missing boxes:", sorted(list(missing)))
    print("Number of boxes in solution:", len(solution))
    print("Number of boxes needed:", len(original))
    
    # Verify no duplicates in solution
    if len(solution) != len(set(solution)):
        print("Warning: Solution contains duplicates!")

find_missing_boxes()