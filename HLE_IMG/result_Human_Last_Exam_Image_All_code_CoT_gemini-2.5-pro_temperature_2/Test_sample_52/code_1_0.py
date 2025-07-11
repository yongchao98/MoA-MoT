import sympy

def calculate_periods():
    """
    Calculates the periods for the four shapes based on geometric analysis.
    The analysis concluded the following periods:
    1) Shape(13, 31, 23): Period = 3
    2) Shape(10, 4, 23, 31): Period = 2
    3) Shape(5, 15, 17, 19, 21, 7): Period = 1
    4) Shape(4, 5, ... 13): Period = 3
    
    This function will print these results.
    """
    
    # Results from the step-by-step analysis
    period1 = 3
    period2 = 2
    period3 = 1
    period4 = 3
    
    # Print the sequence as requested
    print(f"{period1},{period2},{period3},{period4}")

calculate_periods()

# Final Answer must be extracted
final_answer = "3,2,1,3"
# print(f'<<<{final_answer}>>>')