def show_steps():
    val = 98  # After the loop
    print(f"Initial: {val} ({chr(val)})")
    
    val += 10
    print(f"After +10: {val} ({chr(val)})")
    
    val -= 11
    print(f"After -11: {val} ({chr(val)})")
    
    val += 3
    print(f"After +3: {val} ({chr(val)})")
    
    val += 1
    print(f"After +1: {val} ({chr(val)})")

show_steps()