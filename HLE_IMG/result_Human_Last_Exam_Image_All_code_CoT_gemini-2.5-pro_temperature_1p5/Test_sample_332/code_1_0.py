def solve_riddle():
    """
    This script solves the riddle by providing the value for 'Y'.

    The riddle involves these steps:
    1. The image shows Sagan Crater on Mars, named after Carl Sagan.
    2. The riddle links the crater's name to an etymological origin meaning "X of Y".
    3. The surname "Sagan" derives from the Hebrew/Aramaic title "səḡan", meaning deputy or assistant.
    4. This title was notably used as "Səḡan hakKōhănîm", which translates to "Deputy of the Priests".
    5. In the phrase "X of Y" ("Deputy of the Priests"), 'Y' is "Priests".
    """
    
    X = "Deputy"
    Y = "Priests"
    
    print(f"The etymological origin of the name 'Sagan' can be traced to the phrase 'Deputy of the Priests'.")
    print(f"This phrase matches the format 'X of Y'.")
    print(f"In this context, X = '{X}' and Y = '{Y}'.")
    print(f"Therefore, the answer is:")
    print(Y)

solve_riddle()