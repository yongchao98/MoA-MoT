def print_controller_parametrization():
    """
    Prints the parametrized form of the stabilizing controller H_2(s).
    """
    
    nominator = "(s^2 - 1) * K(s) - 1"
    denominator = "s * (1 - K(s))"
    
    print("The set of all proper stabilizing controllers H_2(s) is given by:")
    print("")
    print(f"H_2(s) = {nominator} / {denominator}")
    print("")
    print("where K(s) is any stable and proper rational function.")

if __name__ == '__main__':
    print_controller_parametrization()
