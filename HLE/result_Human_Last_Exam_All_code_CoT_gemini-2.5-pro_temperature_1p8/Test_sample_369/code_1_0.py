def scan_hexameter_line():
    """
    Prints the scansion of the given Latin hexameter line.
    The line is "verum ubi equi atque hominis casu convenit imago."
    D = Dactyl (long, short, short)
    S = Spondee (long, long)
    """
    # The established scansion for this line of Lucretius
    foot1 = "D"
    foot2 = "S"
    foot3 = "D"
    foot4 = "D"
    foot5 = "D"
    foot6 = "S"

    print(f"The scansion of the line is:")
    print(f"{foot1} {foot2} {foot3} {foot4} {foot5} {foot6}")

scan_hexameter_line()