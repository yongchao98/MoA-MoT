def check_missing_boxes():
    original = [48, 293, 96, 382, 284, 61, 122, 234, 211, 109, 202, 133, 179, 293, 182, 168, 313, 101, 43, 262, 100, 231, 159, 280]
    lifted = [382, 293, 313, 293, 48, 284, 280, 43, 262, 234, 109, 231, 211, 101, 202, 182, 168, 179, 159, 133, 122, 100, 96, 61]
    
    # Sort both lists for easier comparison
    original.sort()
    lifted.sort()
    
    # Find missing boxes
    missing = []
    for box in original:
        if lifted.count(box) < original.count(box):
            missing.extend([box] * (original.count(box) - lifted.count(box)))
    
    print("Missing boxes:", missing)
    return missing

print(check_missing_boxes())