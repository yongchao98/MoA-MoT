import sys

def solve():
    """
    This function calculates the minimum number of moves to reverse a sequence of 100 elements
    based on the allowed operations.
    """

    N = 100
    NUM_CLASSES = 5
    ELEMENTS_PER_CLASS = N // NUM_CLASSES

    # Step 1: Explain the concept of congruence classes
    print("Step 1: Understanding the Problem")
    print("The non-adjacent swap (swapping elements with 4 elements between them) is free.")
    print("This operation allows moving an element by 5 positions at no cost.")
    print("This effectively groups the 100 positions into 5 'congruence classes' based on their position modulo 5.")
    print(f"We can reorder all elements within the same class for free.\n")

    # Step 2: Determine where elements need to go
    print("Step 2: Calculating Required Class Shifts")
    print("To reverse the sequence, an element at position 'p' must move to position 'p_final = 101 - p'.")
    print("Let's define classes by 'position % 5'. For multiples of 5, we'll call it Class 5.")
    print("The required class changes are:")
    print("  - Class 1 (p=1,6,..) -> Class 5 (p_final=100,95,..): Shift of -1 (or +4)")
    print("  - Class 2 (p=2,7,..) -> Class 4 (p_final=99,94,..): Shift of +2 (or -3)")
    print("  - Class 3 (p=3,8,..) -> Class 3 (p_final=98,93,..): Shift of  0")
    print("  - Class 4 (p=4,9,..) -> Class 2 (p_final=97,92,..): Shift of -2 (or +3)")
    print("  - Class 5 (p=5,10,..) -> Class 1 (p_final=96,91,..): Shift of +1 (or -4)\n")

    # Step 3: Calculate total shifts needed
    print("Step 3: Calculating Total Work")
    print("A 'move' (adjacent swap) shifts one element's class by +1 and another's by -1.")
    print("The total number of moves is the total number of positive or negative shifts needed.")
    print(f"There are {ELEMENTS_PER_CLASS} elements in each class.")

    print("\nCalculating total positive shifts (work required):")
    # Elements from Class 5 need a +1 shift to get to Class 1
    shift_class_5 = 1
    pos_work_1 = ELEMENTS_PER_CLASS * shift_class_5
    print(f"  - {ELEMENTS_PER_CLASS} elements from Class 5 need a shortest shift of +{shift_class_5}.")

    # Elements from Class 2 need a +2 shift to get to Class 4
    shift_class_2 = 2
    pos_work_2 = ELEMENTS_PER_CLASS * shift_class_2
    print(f"  - {ELEMENTS_PER_CLASS} elements from Class 2 need a shortest shift of +{shift_class_2}.")
    
    total_positive_work = pos_work_1 + pos_work_2
    print(f"\nTotal positive shifts required is the sum of work for these classes.")

    # Final Calculation
    print("\nStep 4: The Final Equation")
    print("The minimum number of moves equals the total positive shifts required.")
    
    # Using 'print' to output each number in the final equation
    sys.stdout.write("(")
    sys.stdout.write(str(ELEMENTS_PER_CLASS))
    sys.stdout.write(" * ")
    sys.stdout.write(str(shift_class_5))
    sys.stdout.write(") + (")
    sys.stdout.write(str(ELEMENTS_PER_CLASS))
    sys.stdout.write(" * ")
    sys.stdout.write(str(shift_class_2))
    sys.stdout.write(") = ")
    sys.stdout.write(str(pos_work_1))
    sys.stdout.write(" + ")
    sys.stdout.write(str(pos_work_2))
    sys.stdout.write(" = ")
    sys.stdout.write(str(total_positive_work))
    sys.stdout.write("\n")
    
    global final_answer
    final_answer = total_positive_work


final_answer = None
solve()
print(f"\n<<<60>>>")