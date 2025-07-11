import sys
import io

# Backup stdout
stdout_backup = sys.stdout
# Redirect stdout to a string buffer
sys.stdout = io.StringIO()

# --- Your code starts here ---

# Given parameters for Snowball Inc.
# You can change these values to test with different scenarios.

# q: number of outstanding shares
q = 1000000

# E: total market value of equity
E = 45000000

# d: total dividends to be distributed in year 1 under the old policy
d = 1500000

# g: annual dividend growth rate
g = 0.05

# --- Calculation ---

# The formula to calculate the per-share ex-dividend price in year 1 (p1) is:
# p1 = (E * (1 + g) - d) / q

p1 = (E * (1 + g) - d) / q

# --- Output ---

# Print the final equation with the numbers plugged in
print(f"The calculation for the per-share ex-dividend price (p1) is:")
print(f"p1 = (E * (1 + g) - d) / q")
print(f"p1 = ({E} * (1 + {g}) - {d}) / {q}")
print(f"p1 = ({E * (1 + g)} - {d}) / {q}")
print(f"p1 = {E * (1 + g) - d} / {q}")
print(f"p1 = {p1}")

# --- Your code ends here ---

# Restore stdout and get the captured output
captured_output = sys.stdout.getvalue()
sys.stdout = stdout_backup

# You can now use the captured_output variable
# For example, print it to the console
print(captured_output)