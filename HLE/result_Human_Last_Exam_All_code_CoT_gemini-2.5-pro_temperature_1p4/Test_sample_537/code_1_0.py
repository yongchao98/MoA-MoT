def solve_equation():
    """
    This function solves the puzzle by constructing and printing the hidden equation.
    The scrambled word "anhaeketn" is an anagram for "TAKE A HEN",
    which translates to the mathematical equation 10 - 1 = 9.
    """
    
    # The numbers derived from the decoded phrase "TAKE A HEN"
    num1 = 10  # from "HEN", which sounds like "TEN"
    num2 = 1   # from "A"
    
    # The operation from the decoded phrase
    operation = "TAKE" # meaning minus
    
    # Calculate the result
    result = num1 - num2
    
    # Print the equation as requested, including each number.
    print(f"The decoded equation is: {num1} - {num2} = {result}")
    print("The numbers in the equation are:")
    print(num1)
    print(num2)
    print(result)

solve_equation()