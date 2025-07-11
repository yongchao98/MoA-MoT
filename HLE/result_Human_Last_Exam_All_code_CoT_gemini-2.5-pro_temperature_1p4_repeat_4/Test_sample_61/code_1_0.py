import sys
from io import StringIO

# This script is designed to print the solution to a legal interpretation question.
# The user wants to know which sections of the UK's Employment Rights Act 1996
# are relevant for determining the procedure for a statutory instrument
# used to adjust section 43FA(3)(c).

# 1. Section 43FA(3)(c): This is the source of the power. It allows the Secretary of
#    State to make regulations to specify who is a "relevant person".
# 2. Section 236(1): This specifies that any power to make regulations under the Act
#    must be exercised by statutory instrument.
# 3. Section 236(2): This establishes the default parliamentary procedure for such
#    statutory instruments, which is the negative resolution procedure.
# 4. Section 236(5): This section is relevant because it lists the exceptions to
#    the default rule in 236(2). One must check this list to confirm that the
#    power under 43FA(3)(c) is not subject to the affirmative procedure. Since
#    43FA(3)(c) is not listed here, the default procedure from 236(2) applies.

# The final list of relevant sections is therefore 43FA(3)(c), 236(1), 236(2), and 236(5).

relevant_sections = "43FA(3)(c),236(1),236(2),236(5)"
print(relevant_sections)
