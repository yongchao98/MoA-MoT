import sys
from io import StringIO

def solve():
  """
  This script identifies the single line of code with the error that causes the
  test case to fail and determines the number of edits to fix it.

  The error is in the loop termination logic. The check `if (feof(stdin)) break;`
  on line 11 occurs after the sums have already been corrupted by an extra loop
  iteration that runs on stale data after `scanf` fails on End-Of-File.

  The most direct error is this faulty check. Removing this line (a single edit
  operation) can be considered a minimal fix that removes the incorrect behavior,
  even though a more complete fix would require restructuring the loop.
  """

  # The faulty line number
  line_number = 11

  # The number of edits to fix it (in this case, deleting the line)
  edit_count = 1

  print(f"{line_number}:{edit_count}")

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = StringIO()

solve()

# Get the captured output
output = captured_output.getvalue().strip()

# Restore stdout
sys.stdout = old_stdout

# Final answer format
# We are asked not to print anything else but the final answer block.
# print(f'The identified error is on line {line_number}, and it can be fixed with {edit_count} edit(s).')
# print('Final Answer:')
print(f'<<<11:1>>>')