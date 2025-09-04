import math

def check_answer():
    """
    Checks the correctness of the provided answer for the light speed in a moving medium problem.
    """

    # The answer from the LLM corresponds to option B.
    # Let's define a function for this formula.
    def formula_B(n, v):
        """Calculates speed based on the formula: (1 + n*v) / (n + v)"""
        # Physical constraints: n >= 1, 0 <= v < 1. So n+v is always positive.
        if n + v == 0:
            return float('inf')
        return (1 + n * v) / (n + v)

    # Test cases with expected results calculated from the fundamental
    # relativistic velocity-addition formula: u = (u' + v) / (1 + u'*v), where u' = 1/n.
    test_cases = [
        {
            "name": "Vacuum Case (n=1)",
            "inputs": {"n": 1.0, "v": 0.5},
            # Expected: In vacuum, speed is always c=1.
            # Formula: (1/1 + 0.5) / (1 + 0.5/1) = 1.0
            "expected": 1.0
        },
        {
            "name": "Stationary Medium Case (v=0)",
            "inputs": {"n": 1.5, "v": 0.0},
            # Expected: In a stationary medium, speed is c/n = 1/n.
            # Formula: (1/1.5 + 0) / (1 + 0) = 1/1.5
            "expected": 1.0 / 1.5
        },
        {
            "name": "Typical Case 1",
            "inputs": {"n": 1.5, "v": 0.5},
            # Expected: u' = 1/1.5; u = (1/1.5 + 0.5) / (1 + (1/1.5)*0.5) = 0.875
            "expected": 0.875
        },
        {
            "name": "Typical Case 2",
            "inputs": {"n": 2.0, "v": 0.1},
            # Expected: u' = 1/2; u = (0.5 + 0.1) / (1 + 0.5*0.1) = 0.6 / 1.05 = 4/7
            "expected": 4.0 / 7.0
        }
    ]

    # Check if the given formula passes all test cases.
    for test in test_cases:
        n = test["inputs"]["n"]
        v = test["inputs"]["v"]
        
        # Check physical constraints of the inputs
        if not (n >= 1):
            return f"Invalid test case: Refractive index n={n} must be >= 1."
        if not (0 <= v < 1):
            return f"Invalid test case: Velocity v={v} must be between 0 and 1 (c)."

        # Calculate result using the formula from the answer
        result = formula_B(n, v)
        expected = test["expected"]

        # Compare with the expected result
        if not math.isclose(result, expected, rel_tol=1e-9):
            reason = f"The answer 'B' is incorrect.\n"
            reason += f"It failed the test '{test['name']}'.\n"
            reason += f"For inputs n={n}, v={v}, the expected result based on the relativistic velocity addition formula is {expected:.6f}.\n"
            reason += f"However, the formula for option 'B' gives {result:.6f}."
            return reason

    # If all tests pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
print(check_answer())