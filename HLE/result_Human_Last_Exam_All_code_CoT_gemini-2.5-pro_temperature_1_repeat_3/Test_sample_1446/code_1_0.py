import collections

def get_critical_exponent_nu():
    """
    Calculates and explains the critical exponent nu for the G4-theory (Ising universality class)
    in different spatial dimensions (d).

    The G4-theoretical framework refers to the standard phi-4 theory used to describe
    continuous phase transitions. The critical exponent nu, which describes the scaling
    of the correlation length (xi ~ |T - T_c|^-nu), is not a single value but
    depends critically on the spatial dimension 'd'.
    """

    print("The value of the critical exponent ν for the φ⁴-theory (Ising universality class) depends on the spatial dimension d:")
    print("-" * 80)

    # Using an ordered dictionary to present the dimensions in a logical sequence
    nu_values = collections.OrderedDict()
    
    # For d >= 4 (upper critical dimension), mean-field theory is exact.
    nu_values["d >= 4"] = {
        "value": 0.5,
        "type": "Exact (Mean-Field Theory)",
        "description": "This is the upper critical dimension. Fluctuations are irrelevant, and mean-field theory holds."
    }
    
    # For d = 3, determined by high-precision numerical and analytical methods.
    nu_values["d = 3"] = {
        "value": 0.629971,
        "type": "High-precision result (Renormalization Group, Bootstrap)",
        "description": "This is the most physically relevant non-trivial case, describing 3D Ising-like systems."
    }

    # For d = 2, from the exact solution of the 2D Ising model.
    nu_values["d = 2"] = {
        "value": 1.0,
        "type": "Exact (Onsager Solution)",
        "description": "The 2D case is exactly solvable and serves as a benchmark for statistical mechanics."
    }

    # Print the information for each dimension.
    # The prompt requests to "output each number in the final equation".
    # We interpret this as clearly stating the value for each case.
    for dim, data in nu_values.items():
        print(f"For dimension {dim}:")
        # In finite-size scaling, the correlation length ξ scales with system size L.
        # Near the critical temperature T_c, ξ also scales with the reduced temperature t = |T-T_c|/T_c.
        # The relationship is ξ ~ |t|^(-ν).
        print(f"  The scaling equation involves ν = {data['value']}")
        print(f"  This value is: {data['type']}")
        print(f"  Description: {data['description']}\n")

if __name__ == "__main__":
    get_critical_exponent_nu()
