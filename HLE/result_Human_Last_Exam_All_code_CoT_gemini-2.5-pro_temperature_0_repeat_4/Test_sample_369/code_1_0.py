def scan_hexameter_line():
    """
    Prints the scansion of the specified Latin hexameter line.
    D = Dactyl (— U U)
    S = Spondee (— —)
    """
    line = "verum ubi equi atque hominis casu convenit imago."
    foot_1 = "D"
    foot_2 = "S"
    foot_3 = "D"
    foot_4 = "S"
    foot_5 = "D"
    foot_6 = "S"
    
    scansion = f"{foot_1} {foot_2} {foot_3} {foot_4} {foot_5} {foot_6}"
    
    print(f"The line is: \"{line}\"")
    print(f"The scansion is: {scansion}")

scan_hexameter_line()