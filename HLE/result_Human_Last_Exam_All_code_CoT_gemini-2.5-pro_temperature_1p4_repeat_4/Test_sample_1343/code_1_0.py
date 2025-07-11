def solve():
    """
    This script simulates the high-level logic of the provided RASPy code.
    The RASPy function 'v' checks for a special condition using function 'u'.
    Function 'u' counts the occurrences of the digit '7' in three specific
    3-character segments of the input string: sop[0:3], sop[3:6], and the last 3 characters.
    If the total count of '7's is exactly 1, a fixed string is output.
    Otherwise, the string is treated as an addition problem (e.g., "num1+num2"),
    and the sum is calculated and returned.
    """
    inputs = ["734107+4295754", "5429141+142196"]
    outputs = []

    for sop in inputs:
        # Simulate function u(sop) to count the '7's
        n = len(sop)
        chunk1 = sop[0:3]
        chunk2 = sop[3:6]
        chunk3 = sop[n-3:n]
        
        sevens_count = chunk1.count('7') + chunk2.count('7') + chunk3.count('7')
        
        # Simulate the main logic of function v(sop)
        if sevens_count == 1:
            # This branch is not taken by the given inputs
            result_str = "get pwned!..."
            outputs.append(result_str)
            print(f"Input '{sop}' has 1 seven in key segments. Magic output: {result_str}")
        else:
            # Perform addition as the '7's count is not 1
            try:
                parts = sop.split('+')
                num1 = int(parts[0])
                num2 = int(parts[1])
                result_val = num1 + num2
                outputs.append(str(result_val))
                # Print each number in the final equation as requested
                print(f"{num1} + {num2} = {result_val}")
            except (ValueError, IndexError):
                # Handle cases where the input is not a valid addition expression
                outputs.append("Error: Invalid input format")

solve()