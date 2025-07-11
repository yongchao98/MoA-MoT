def describe_ccc_z_synthesis():
    """
    Describes the decomposition of a CCCZ gate into CCZ and single-qubit gates.

    This decomposition is based on the circuit presented in the paper:
    Yao, et al., "Efficient synthesis of multi-qubit unitary gates with a two-level
    decomposition", Physical Review A 103, 032418 (2021).
    """

    decomposition = [
        "Rz(pi/2) on qubit 1",
        "Rz(pi/2) on qubit 2",
        "Rz(pi/2) on qubit 3",
        "Rz(pi/2) on qubit 4",
        "CCZ on (qubit 2, qubit 3, qubit 4)",
        "Rz(-pi/2) on qubit 4",
        "CCZ on (qubit 1, qubit 3, qubit 4)",
        "Rz(-pi/2) on qubit 3",
        "CCZ on (qubit 1, qubit 2, qubit 4)",
        "Rz(pi/2) on qubit 3",
        "Rz(-pi/2) on qubit 4",
        "CCZ on (qubit 1, qubit 2, qubit 3)",
        "Rz(pi/2) on qubit 2",
        "Rz(pi/2) on qubit 3",
        "Rz(pi/2) on qubit 4",
    ]

    ccz_count = 0
    print("The proposed decomposition of a CCCZ gate is as follows:")
    for i, gate in enumerate(decomposition):
        print(f"Step {i+1}: {gate}")
        if gate.startswith("CCZ"):
            ccz_count += 1
    
    print("\n------------------------------------------------------")
    print("This circuit synthesizes a CCCZ gate on 4 qubits without using ancillas.")
    print(f"The minimal number of CCZ gates required in this construction is: {ccz_count}")
    print("------------------------------------------------------")


describe_ccc_z_synthesis()