import sys
import io

# Redirect stdout to a string buffer to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve():
    """
    This function provides the solution based on the visual analysis of the 16 plots.
    Each character corresponds to a plot from 1 to 16.
    - B/b: Blue vortex is unique (twice/half strength)
    - G/g: Green vortex is unique (twice/half strength)
    - R/r: Red vortex is unique (twice/half strength)
    """
    
    # Analysis results for each plot from 1 to 16
    # Plot 1: Blue is weak (b)
    # Plot 2: Red is strong (R)
    # Plot 3: Green is strong (G)
    # Plot 4: Blue is strong (B)
    # Plot 5: Blue is strong (B)
    # Plot 6: Green is weak (g)
    # Plot 7: Red is strong (R)
    # Plot 8: Red is weak (r)
    # Plot 9: Green is weak (g)
    # Plot 10: Red is strong (R)
    # Plot 11: Blue is strong (B)
    # Plot 12: Blue is weak (b)
    # Plot 13: Green is weak (g)
    # Plot 14: Blue is weak (b)
    # Plot 15: Red is strong (R)
    # Plot 16: Green is weak (g)
    
    result = "bRGBBgRrGBbgRbRg"
    print(result)

solve()

# Restore stdout and get the captured output
sys.stdout = old_stdout
output = captured_output.getvalue()

# Print the final result in the required format
print(output)
print(f'<<<{output.strip()}>>>')